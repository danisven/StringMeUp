#!/usr/bin/env python3

__version__ = "0.2.0"  # Major bug fix in confidence calculation (ambiguous k-mers)

import argparse
import operator
import logging
import gzip
import sys
from stringmeup import taxonomy
from dataclasses import dataclass
from os import path
from typing import Dict, Set, Tuple, Optional, List, TextIO

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d [%H:%M:%S]')
log = logging.getLogger(path.basename(__file__))

# TODO: make sure confidence_threshold is between 0 and 1
# TODO: For the verbose output, also output (1) the number of kmers that hit in total, (2) the number of non-ambiguous kmers (queried).


@dataclass
class ReadClassification:
    current_node: int = None
    original_conf: float = None
    recalculated_conf: float = None
    original_taxid: int = None
    reclassified_taxid: int = None
    original_rank_code: str = None
    reclassified_rank_code: str = None
    original_name: str = None
    reclassified_name: str = None
    reclassified_distance: int = None
    id: str = None
    length: str = None
    kmer_string: str = None
    classified: bool = False
    max_confidence = None
    minimizer_hit_groups = None
    ambiguous_kmers: int = 0  # Count of ambiguous k-mers (A:X in Kraken2 output)


@dataclass
class ReportNode:
    ratio: str
    hits_at_clade: int
    hits_at_node: int
    rank_code: str
    rank_depth: int
    node_taxid: int
    name: str
    offset: int


def validate_input_file(putative_classifications_file: str, verbose_input: bool, minimum_hit_groups: Optional[int], paired_input: bool) -> None:
    """
    Perform simple validation of the input file.
    """

    log.debug('Validating input classifications file.')

    if not path.isfile(putative_classifications_file):
        log.error('Cannot find the specified file ({file}).'.format(
            file=putative_classifications_file))
        sys.exit()

    with read_file(putative_classifications_file) as f:
        line = f.readline()
        line_proc = line.strip()
        line_proc = line_proc.split('\t')

        # The following should be the case of a Kraken 2 output file
        # First, check so the number of columns in the input file conforms to the expected number
        if not verbose_input:
            num_cols = len(line_proc) == 5  # original type of kraken2 output file
        else:
            num_cols = len(line_proc) == 6  # 6 columns if the output was produced with the verbose version of kraken2 that outputs minimizer hit groups

        # Line must start with C or U (as in Classified/unclassified)
        line_start = line_proc[0] in ['U', 'C']

        # If the data is paired
        if paired_input:
            # Length field (index 3) should still be an integer
            try:
                int(line_proc[3])
            except (ValueError, TypeError, IndexError):
                data_col_1 = False
            else:
                data_col_1 = True

            # K-mer string (last column) must contain |:| delimiter with data on both sides
            if "|" in line_proc[-1]:
                data_col_2 = len(line_proc[-1].split('|:|')) == 2
            else:
                data_col_2 = False

        # If the input is from single end reads, atleast the read length column (3rd) must be an int
        else:
            try:
                int(line_proc[3])
            except (ValueError, TypeError, IndexError):
                data_col_1 = False
            else:
                data_col_1 = True

            # And the last column should contain colons between kmer/taxon pairs
            if ":" in line_proc[4]:
                data_col_2 = True
            else:
                data_col_2 = False

        if num_cols and line_start and data_col_1 and data_col_2:
            log.debug('Validation OK.')
            return
        else:
            log.error('The classifications file is malformatted.')
            log.debug('First line of input: {}'.format(line))
            log.debug('num_cols: {}'.format(num_cols))
            log.debug('line_start: {}'.format(line_start))
            log.debug('data_col_1: {}'.format(data_col_1))
            log.debug('data_col_1: {}'.format(data_col_2))
            sys.exit()


def is_paired_input(classifications_file: str) -> bool:
    """
    Returns true if input file appears to contain paired read data.
    """
    with read_file(classifications_file) as f:
        line = f.readline()
        line_proc = line.strip()
        line_proc = line_proc.split('\t')

        # If k-mer string (index 4) contains a pipe character "|", the data is paired
        if "|" in line_proc[4]:
            return True


def is_verbose_input(classifications_file: str) -> bool:
    """
    Returns true if input file consists of 6 columns instead of 5.
    """

    with read_file(classifications_file) as f:
        line = f.readline()
        line_proc = line.strip()
        line_proc = line_proc.split('\t')
        if len(line_proc) == 6:
            return True
        else:
            return False


def process_kmer_string(kmer_info_string: str, paired_input: bool) -> Tuple[Dict[int, int], int]:
    """
    Process a kmer info string (last column of a Kraken 2 output file), so that
    we get a dictionary mapping of tax_ids to total sum of kmer hits, plus the count
    of ambiguous k-mers.

    Returns:
        Tuple of (taxa_kmer_dict, ambiguous_kmer_count)
        - taxa_kmer_dict: {tax_id_#1: X kmer hits, tax_id_#2: Y kmer hits, ...}
        - ambiguous_kmer_count: Total count of ambiguous k-mers (marked with 'A' in output)

    CRITICAL: Ambiguous k-mers MUST be counted separately and added to the total denominator
    to match Kraken2's confidence calculation. Kraken2 includes ambiguous k-mers in total_kmers
    even though they cannot be queried against the database.
    """
    kmer_info_string = kmer_info_string.split()

    # Kraken2 classifications file for paired data contain the "|:|" delimiter
    if paired_input:
        kmer_info_string.remove('|:|')

    # Process all "taxa":"num_kmer" string pairs, including ambiguous k-mers marked with 'A'
    # Structural delimiters like '-:-' (reading frame borders) are skipped
    taxa_kmer_dict = {}
    ambiguous_kmer_count = 0

    for kmer_info in kmer_info_string:
        parts = kmer_info.split(':')

        # Skip malformed entries
        if len(parts) != 2:
            continue

        tax_id_str, count_str = parts

        # Skip reading frame borders (used in translated search)
        if tax_id_str == '-':
            continue

        try:
            count = int(count_str)
        except ValueError:
            continue

        # Handle ambiguous k-mers separately
        if tax_id_str == 'A':
            ambiguous_kmer_count += count
        else:
            # Regular tax_id
            try:
                tax_id = int(tax_id_str)
                if tax_id not in taxa_kmer_dict:
                    taxa_kmer_dict[tax_id] = count
                else:
                    taxa_kmer_dict[tax_id] += count
            except ValueError:
                # Skip if tax_id is not a valid integer
                continue

    return taxa_kmer_dict, ambiguous_kmer_count


def _calculate_max_confidence(taxa_kmer_dict: Dict[int, int], total_kmer_hits: int) -> float:
    """
    Calculate the maximum possible confidence for a read.
    This is the sum of all k-mers assigned to taxa (excluding 0) divided by total k-mers.
    """
    assigned_taxa_set = set([tax_id for tax_id in taxa_kmer_dict.keys() if tax_id != 0])
    total_hits = sum(taxa_kmer_dict[tax_id] for tax_id in assigned_taxa_set)
    return total_hits / total_kmer_hits


def _check_early_termination(max_confidence: float, confidence_threshold: float,
                              verbose_input: bool, minimizer_hit_groups: int,
                              minimum_hit_groups: int) -> bool:
    """
    Check if the read is "doomed to fail" classification based on:
    1. Maximum possible confidence is below threshold
    2. Minimizer hit groups (if using verbose input) is below minimum

    Returns True if the read cannot possibly meet classification criteria.
    """
    if max_confidence < confidence_threshold:
        return True

    if verbose_input and minimizer_hit_groups < minimum_hit_groups:
        return True

    return False


def _get_cached_lineage(tax_id: int, taxa_lineages: Dict[int, List[int]],
                        taxonomy_tree: taxonomy.TaxonomyTree) -> List[int]:
    """
    Get the lineage for a tax_id, using cache if available, otherwise fetch and cache it.
    """
    if tax_id in taxa_lineages:
        return taxa_lineages[tax_id]
    else:
        lineage = taxonomy_tree.get_lineage([tax_id])[tax_id]
        # Save lineage so we don't have to get it from taxonomy_tree more than once.
        # Keep as list to maintain consistency with taxonomy_tree.get_lineage() return type.
        taxa_lineages[tax_id] = lineage
        return lineage


def _calculate_clade_confidence(current_node: int, assigned_taxa_set: Set[int],
                                 taxa_kmer_dict: Dict[int, int], total_kmer_hits: int,
                                 taxa_lineages: Dict[int, List[int]],
                                 taxonomy_tree: taxonomy.TaxonomyTree) -> Tuple[float, Set[int], Set[int]]:
    """
    Calculate confidence at current node by summing k-mers from all descendants.

    Returns:
        - confidence score at this node
        - set of tax_ids found in this clade
        - updated assigned_taxa_set with clade taxa removed
    """
    num_hits_within_clade = 0
    taxa_in_clade = set()

    for tax_id in assigned_taxa_set:
        lineage = _get_cached_lineage(tax_id, taxa_lineages, taxonomy_tree)

        # If current_node is in the lineage of tax_id, then tax_id is a descendant
        if current_node in lineage:
            taxa_in_clade.add(tax_id)
            num_hits_within_clade += taxa_kmer_dict[tax_id]

    # Remove counted taxa from the search set
    updated_assigned_taxa_set = assigned_taxa_set - taxa_in_clade if taxa_in_clade else assigned_taxa_set

    # Early exit optimization: if all taxa counted, sum all k-mers
    if not updated_assigned_taxa_set:
        num_hits_within_clade = sum(taxa_kmer_dict[tax_id] for tax_id in taxa_kmer_dict.keys() if tax_id != 0)

    confidence = num_hits_within_clade / total_kmer_hits
    return confidence, taxa_in_clade, updated_assigned_taxa_set


def reclassify_read(read: ReadClassification, confidence_threshold: float, taxonomy_tree: taxonomy.TaxonomyTree, verbose_input: bool, minimum_hit_groups: int, taxa_lineages: Dict[int, List[int]], paired_input: bool) -> Tuple[ReadClassification, Dict[int, List[int]]]:
    """
    Reclassifies a read by traversing up the taxonomy tree until confidence threshold is met.

    Algorithm:
    1. Calculate confidence = (k-mers in clade) / (total k-mers including ambiguous)
    2. If confidence >= threshold, classify to current node
    3. Otherwise, move to parent node and repeat
    4. Stop at root or when confidence threshold is met

    CRITICAL: Must match Kraken2's confidence calculation, which includes ambiguous k-mers
    in the denominator even though they cannot be queried against the database.
    """
    # Process the kmer string into a dict of {tax_id: #kmers} key, value pairs
    # CRITICAL: ambiguous_kmer_count must be added to total for correct confidence calculation
    taxa_kmer_dict, ambiguous_kmer_count = process_kmer_string(read.kmer_string, paired_input)

    # Store ambiguous count in read object
    read.ambiguous_kmers = ambiguous_kmer_count

    # Make the current node the same as the original classification
    read.current_node = read.original_taxid

    # Total number of k-mers - MUST include ambiguous to match Kraken2
    # This represents ALL k-mers in the read, not just those queried
    total_kmer_hits = sum(taxa_kmer_dict.values()) + ambiguous_kmer_count

    # Handle edge case where read has no k-mers (empty or all ambiguous)
    if total_kmer_hits == 0:
        read.reclassified_taxid = 0
        read.classified = False
        read.original_conf = 0.0
        read.recalculated_conf = 0.0
        read.max_confidence = 0.0
        return read, taxa_lineages

    # Only interested in tax_ids that are in the database (exclude 0)
    assigned_taxa_set = set([tax_id for tax_id in taxa_kmer_dict.keys() if tax_id != 0])

    # Check if classification is possible (early termination optimization)
    max_confidence = _calculate_max_confidence(taxa_kmer_dict, total_kmer_hits)
    read.max_confidence = max_confidence

    doomed_to_fail = _check_early_termination(
        max_confidence, confidence_threshold, verbose_input,
        read.minimizer_hit_groups or 0, minimum_hit_groups
    )

    # Traverse up the taxonomy tree until confidence threshold is met
    while not read.classified:
        # Calculate confidence at current node
        confidence, taxa_in_clade, assigned_taxa_set = _calculate_clade_confidence(
            read.current_node, assigned_taxa_set, taxa_kmer_dict,
            total_kmer_hits, taxa_lineages, taxonomy_tree
        )

        read.recalculated_conf = confidence

        # Set original confidence on first iteration
        if not read.original_conf:
            read.original_conf = read.recalculated_conf

            # Exit early if we know classification will fail
            if doomed_to_fail:
                read.recalculated_conf = max_confidence
                read.reclassified_taxid = 0
                return read, taxa_lineages

        # Check if confidence threshold is met
        if read.recalculated_conf >= confidence_threshold:
            read.classified = True
            read.reclassified_taxid = read.current_node
            return read, taxa_lineages

        # Check if we've reached the root (can't go higher)
        elif read.current_node == 1:
            read.classified = False
            read.reclassified_taxid = 0
            return read, taxa_lineages

        # Move to parent node and continue
        else:
            read.current_node = taxonomy_tree.get_parent([read.current_node])[read.current_node]


def get_kraken2_report_content(tax_reads: Dict[str, Dict[int, int]], taxonomy_tree: taxonomy.TaxonomyTree, total_reads: int) -> List[ReportNode]:
    """
    First calculates the cumulative clade read count (i.e. for each node, how
    many reads are classified to the clades rooted at that node).

    Then, sorts the nodes in the order they will be printed in the report file.
    The sorting works like this: perform a depth first search of the taxonomy,
    starting at the root. At each node, continue with the depth first seach
    in the order of highest to lowest cumulative clade read counts among the
    child nodes.

    There is probably a way of merging the internal functions dfs_ccrc and
    dfs_sort, but I was in a hurry and this is what it is.

    total_reads: total reads in the kraken output file.
    """
    report_node_list = []

    # Depth First Search, cumulative clade read count
    def dfs_ccrc(visited_nodes, taxonomy_tree, node_taxid):
        if node_taxid not in visited_nodes:
            visited_nodes.add(node_taxid)

            # The cumulative clade read count for the clade rooted at current
            # node will start at the number of hits at that node:
            if node_taxid in tax_reads['hits_at_node']:
                hits_at_node = tax_reads['hits_at_node'][node_taxid]
            else:
                hits_at_node = 0
            clade_read_count = hits_at_node

            # Recursively search the children in the clade rooted at the
            # current node:
            for child in taxonomy_tree.get_children([node_taxid])[node_taxid]:
                clade_read_count += dfs_ccrc(visited_nodes, taxonomy_tree, child)

            # Only save the node if its clade read count is !=0:
            if clade_read_count > 0:
                tax_reads['hits_at_clade'][node_taxid] = clade_read_count

            return clade_read_count  # This return ends up in the above for loop

    # Depth First Search, sorting for the hierarchy of the output report.
    # Moves down the tree and adds the nodes with highest cumulative clade
    # read count to report_node_list
    def dfs_sort(visited_nodes, taxonomy_tree, node_taxid, offset):
        if node_taxid not in visited_nodes:
            visited_nodes.add(node_taxid)

            # Find the cumluative read counts for all children of current node:
            children_cumulative_counts = {}
            for child in taxonomy_tree.get_children([node_taxid])[node_taxid]:
                if child in tax_reads['hits_at_clade']:
                    children_cumulative_counts[child] = tax_reads['hits_at_clade'][child]

            # Get the hits to this node (to be output in the report):
            if node_taxid in tax_reads['hits_at_node']:
                hits_at_node = tax_reads['hits_at_node'][node_taxid]
            else:
                hits_at_node = 0

            # Get some information that will go into the report file.
            # 1) total number of reads at this node and at the clade rooted here (column 2 in output):
            hits_at_clade = tax_reads['hits_at_clade'][node_taxid]
            # 2) the ratio (column 1 in output):
            ratio_classified2clade = hits_at_clade / total_reads * 100
            # 3) the rank code and depth (column 4 in output):
            rank_tuple = taxonomy_tree.get_rank_code([node_taxid])[node_taxid]
            # 4) scientific name of the node (column 6 of output):
            name = taxonomy_tree.get_name([node_taxid])[node_taxid]

            # Construct the dataclass instance that holds the information
            # about this node that is printed to the report file:
            report_node = ReportNode(
                ratio="{0:.2f}".format(ratio_classified2clade),
                hits_at_clade=hits_at_clade,
                hits_at_node=hits_at_node,
                rank_code=rank_tuple.rank_code,
                rank_depth=rank_tuple.rank_depth,
                node_taxid=node_taxid,
                offset=offset,
                name=name)

            # Append it to the list. The order of the elements in the list is
            # the order the nodes will be printed.
            report_node_list.append(report_node)

            if len(children_cumulative_counts) > 0:
                # sorted_by_ccrc is a list of tuples [(tax_id#1, ccrc), (tax_id#2, ccrc)...]
                sorted_by_ccrc = sorted(
                    children_cumulative_counts.items(),
                    key=operator.itemgetter(1),
                    reverse=True)

                # We are going one level deeper in the taxonomy, increment
                # the offset:
                offset += 1

                for child_tuple in sorted_by_ccrc:
                    child_taxid = child_tuple[0]
                    dfs_sort(visited_nodes, taxonomy_tree, child_taxid, offset)

    # Depth first search to get cumulative read counts for all clades:
    visited_nodes = set()
    log.info('Calculating cumulative clade read counts...')
    dfs_ccrc(visited_nodes, taxonomy_tree, 1)

    # Depth first search to sort the output according to largest cumulative
    # read count:
    visited_nodes = set()
    log.info('Sorting the order of the output in the report file...')
    dfs_sort(visited_nodes, taxonomy_tree, 1, 0)

    # Make sure to add the unclassified row that goes at the very top of the
    # kraken 2 report:
    num_unclassified_reads = total_reads - tax_reads['hits_at_clade'][1]
    ratio = num_unclassified_reads / total_reads * 100
    unclassified_node = ReportNode(
        ratio="{0:.2f}".format(ratio),
        hits_at_clade=num_unclassified_reads,
        hits_at_node=num_unclassified_reads,
        rank_code='U',
        rank_depth=0,
        node_taxid=0,
        name='unclassified',
        offset=0)

    report_node_list.insert(0, unclassified_node)

    return report_node_list


def format_kraken2_report_row(report_node: ReportNode) -> str:
    """
    Formats the row that will be output in the kraken 2 style report. Input
    is an instance of ReportNode.
    """
    offset = 2 * report_node.offset * ' '
    name = offset + report_node.name

    if report_node.rank_depth == 0:
        rank_depth = ''
    else:
        rank_depth = str(report_node.rank_depth)

    rank_code = report_node.rank_code + rank_depth
    report_row = '\t'.join([
        str(report_node.ratio),
        str(report_node.hits_at_clade),
        str(report_node.hits_at_node),
        rank_code,
        str(report_node.node_taxid),
        name])

    return report_row


def make_kraken2_report(tax_reads: Dict[str, Dict[int, int]], taxonomy_tree: taxonomy.TaxonomyTree, total_reads: int, output_report: Optional[str]) -> None:
    """
    Gets the information that should be printed from
    get_kraken2_report_content. Formats the information and prints it to file
    or stdout.
    """
    log.info('Creating report...')

    # Get the data to print
    report_node_list = get_kraken2_report_content(
        tax_reads, taxonomy_tree, total_reads)

    # If the output should go to file
    if output_report:
        with open(output_report, 'w') as f:
            for node in report_node_list:

                # Format the output
                report_row = format_kraken2_report_row(node)
                f.write(report_row + '\n')

            log.info('Report saved in {}.'.format(output_report))

    # Otherwise, print to stdout
    else:
        for node in report_node_list:

            # Format the output
            report_row = format_kraken2_report_row(node)
            sys.stdout.write(report_row + '\n')


def get_verbose_output(read: ReadClassification, taxonomy_tree: taxonomy.TaxonomyTree) -> ReadClassification:
    """
    Gets more information about the reclassification of a read. This is for
    the output_verbose option.

    read is an instance of ReadClassification.
    """
    # Variable renaming to make things more readable
    new_taxid = read.reclassified_taxid
    old_tax_id = read.original_taxid

    # Stuff that is conditional if we have a classified read or not.
    # TaxonomyTree doesn't cope with tax_id=0
    if read.classified:
        distance = taxonomy_tree.get_distance(new_taxid, old_tax_id)
        new_rank_tuple = taxonomy_tree.get_rank_code([new_taxid])[new_taxid]
        new_rank_depth = new_rank_tuple.rank_depth if new_rank_tuple.rank_depth != 0 else ''
        new_rank_code = str(new_rank_tuple.rank_code) + str(new_rank_depth)
        new_rank_name = taxonomy_tree.get_name([new_taxid])[new_taxid]
    else:
        distance = 'NaN'
        new_rank_code = 'U'
        new_rank_name = 'unclassified'

    # Distance information
    read.reclassified_distance = distance

    # Rank information
    old_rank_tuple = taxonomy_tree.get_rank_code([old_tax_id])[old_tax_id]
    old_rank_depth = old_rank_tuple.rank_depth if old_rank_tuple.rank_depth != 0 else ''
    old_rank_code = str(old_rank_tuple.rank_code) + str(old_rank_depth)
    read.original_rank_code = old_rank_code
    read.reclassified_rank_code = new_rank_code

    # Scientific name information
    old_rank_name = taxonomy_tree.get_name([old_tax_id])[old_tax_id]
    read.original_name = old_rank_name
    read.reclassified_name = new_rank_name

    return read


def create_read(kraken2_read: str, verbose_input: bool = False) -> ReadClassification:
    """
    Creates an instance of ReadClassification dataclass, that holds
    information about the read and its classification.
    """
    # Process the read string so that its elements go into a list
    read_pair_proc = kraken2_read.strip()
    read_pair_proc = read_pair_proc.split('\t')

    # Create the read object
    read = ReadClassification(
        original_taxid=int(read_pair_proc[2]),
        id=read_pair_proc[1],
        length=read_pair_proc[3],
        kmer_string=read_pair_proc[-1])

    if verbose_input:
        read.minimizer_hit_groups = int(read_pair_proc[4])

    return read


def main_loop(f_handle, tax_reads_dict, taxonomy_tree, args, report_frequency, taxa_lineages, paired_input, verbose_input=False, o_handle=None, v_handle=None):
    """
    f_handle: classifications input file to read from.
    o_handle: output_classifications file to write to.
    v_handle: output_verbose file to write to.
    """
    def write_read_output(read):
        # read is an instance of ReadClassification
        classification = 'C' if read.classified else 'U'
        row_items = [
            classification,
            read.id,
            read.reclassified_taxid,
            read.length,
            read.kmer_string]

        if verbose_input:
            row_items.insert(4, read.minimizer_hit_groups)

        row_string = '\t'.join([str(x) for x in row_items]) + '\n'
        _ = o_handle.write(row_string)  # gzip write fnc returns output, therefore send to "_"

    def write_verbose_output(read):
        # read is an instance of ReadClassification
        row_items = [
            read.id,
            read.length,
            read.reclassified_distance,
            read.original_taxid,
            read.reclassified_taxid,
            "{0:.2f}".format(read.original_conf),
            "{0:.2f}".format(read.recalculated_conf),
            "{0:.2f}".format(read.max_confidence),
            read.original_rank_code,
            read.reclassified_rank_code,
            read.original_name,
            read.reclassified_name,
            read.kmer_string]

        if verbose_input:
            row_items.insert(2, read.minimizer_hit_groups)

        row_string = '\t'.join([str(x) for x in row_items]) + '\n'
        _ = v_handle.write(row_string)  # gzip write fnc returns output, therefore send to "_"

    # Parse the input file, read per read
    i = 0
    for read_pair in f_handle:

        # Only working with classified reads:
        if read_pair.startswith('C'):

            # Make an instance of ReadClassification to hold information
            # about the read and its classification
            read = create_read(read_pair, verbose_input)

            # Validate that classified reads have a non-zero tax_id
            if read.original_taxid == 0:
                log.warning('Read {} is marked as classified but has tax_id=0. Skipping.'.format(read.id))
                continue

            # Reclassify the read pair based on confidence
            read, taxa_lineages = reclassify_read(
                read,
                args.confidence_threshold,
                taxonomy_tree,
                verbose_input,
                args.minimum_hit_groups,
                taxa_lineages,
                paired_input)

            # Counter for number of reads per taxon/node
            if read.reclassified_taxid in tax_reads_dict['hits_at_node']:
                tax_reads_dict['hits_at_node'][read.reclassified_taxid] += 1
            else:
                tax_reads_dict['hits_at_node'][read.reclassified_taxid] = 1

            # Write the reclassified reads to file
            if o_handle:
                if read.classified or args.keep_unclassified:
                    write_read_output(read)

            # Write verbose output about the reclassification
            if v_handle:
                read = get_verbose_output(read, taxonomy_tree)
                write_verbose_output(read)

        else:
            # Change here if you want to keep reads from the input file
            # that were initially unclassified.
            pass

        # Keep track of progress
        i += 1
        if i > 0 and i % report_frequency == 0:
            log.info('Processed {} reads...'.format(i))

    log.info('Done processing reads. They were {} in total.'.format(i))

    # Output a report file
    make_kraken2_report(tax_reads_dict, taxonomy_tree, i, args.output_report)  # i is used to calculate the ratio of classified reads (col 1 in output file).


def read_file(filename: str) -> TextIO:
    """
    Wrapper to read either gzipped or ordinary text file input.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


def write_file(filename: str, gz_output: bool) -> TextIO:
    """
    Wrapper to write either gzipped or ordinary text file output.
    """
    if gz_output:
        return gzip.open(filename, 'wt')
    else:
        return open(filename, 'w')


def get_arguments() -> argparse.Namespace:
    """
    Wrapper function to get the command line arguments. Inserting this piece of code
    into its own function for conda compatibility.
    """

    parser = argparse.ArgumentParser(
        prog='StringMeUp',
        usage='stringmeup --names <FILE> --nodes <FILE> [--output_report <FILE>] [--output_classifications <FILE>] [--output_verbose <FILE>] [--keep_unclassified] [--minimum_hit_groups INT] [--gz_output] [--help] confidence classifications',
        description='A post-processing tool to reclassify Kraken 2 output based on the confidence score and/or minimum minimizer hit groups.')
    parser.add_argument(
        'confidence_threshold',
        metavar='confidence',
        type=float,
        help='The confidence score threshold to be used in reclassification [0-1].')
    parser.add_argument(
        'original_classifications_file',
        metavar='classifications',
        type=str,
        help='Path to the Kraken 2 output file containing the individual read classifications.')
    parser.add_argument(
        '--output_report',
        metavar='FILE',
        type=str,
        help='File to save the Kraken 2 report in.')
    parser.add_argument(
        '--output_classifications',
        metavar='FILE',
        type=str,
        help='File to save the Kraken 2 read classifications in.')
    parser.add_argument(
        '--keep_unclassified',
        action='store_true',
        help='Specify if you want to output unclassified reads in addition to classified reads. NOTE(!): This script will always discard reads that are unclassified in the classifications input file, this flag will just make sure to keep previously classified reads even if they are reclassified as unclassified by this script. TIP(!): Always run Kraken2 with no confidence cutoff.')
    parser.add_argument(
        '--output_verbose',
        metavar='FILE',
        type=str,
        help='File to send verbose output to. This file will contain, for each read, (1) original classification, (2) new classification, (3) original confidence, (4), new confidence (5), original taxa name (6), new taxa name, (7) original rank, (8) new rank, (9) distance travelled (how many nodes was it lifted upwards in the taxonomy).')
    parser.add_argument(
        '--names',
        metavar='FILE',
        required=True,
        help='Taxonomy names dump file (names.dmp)')
    parser.add_argument(
        '--nodes',
        metavar='FILE',
        required=True,
        help='Taxonomy nodes dump file (nodes.dmp)')
    parser.add_argument(
        '--minimum_hit_groups',
        metavar='INT',
        type=int,
        help='The minimum number of hit groups a read needs to be classified. NOTE: You need to supply a classifications file (kraken2 output) that contain the "minimizer_hit_groups" column.')
    parser.add_argument(
        '--gz_output',
        action='store_true',
        help='Set this flag to output <output_classifications> and <output_verbose> in gzipped format (will add .gz extension to the filenames).'
    )
    args = parser.parse_args()

    return args


def stringmeup() -> None:

    # Get the CL arguments
    args = get_arguments()

    # Validate confidence threshold is in valid range [0, 1]
    if not 0 <= args.confidence_threshold <= 1:
        log.error('Confidence threshold must be between 0 and 1, got {:.2f}'.format(args.confidence_threshold))
        sys.exit(1)

    # Some initial setup
    taxa_lineages = {}
    report_frequency = 10000000  # Will output progress every nth read
    tax_reads_dict = {'hits_at_node': {}, 'hits_at_clade': {}}

    # Was the input generated with https://github.com/danisven/kraken2 ?
    verbose_input = is_verbose_input(args.original_classifications_file)

    # If so, output warnings if input doesn't contain minimizer hit groups
    if verbose_input:
        log.info("The input file appears to contain a column for minimizer hit groups.")

        if not args.minimum_hit_groups:
            log.warning('You didn\'t specify --minimum_hit_groups.')
            log.warning('Will NOT reclassify based on minimizer hit groups, setting minimum_hit_groups=1 (lowest possible setting).')
            args.minimum_hit_groups = 1

    # If not, output warnings if user was planning to reclassify with minimizer hit groups
    else:
        if args.minimum_hit_groups:
            log.warning('You specified minimum_hit_groups={}, but did not supply an input file with a column for minimizer hit groups.'.format(args.minimum_hit_groups))
            log.warning('Will NOT reclassify based on minimizer hit groups.')
            args.minimum_hit_groups = None

    # Check if the input data is paired or not
    paired_input = is_paired_input(args.original_classifications_file)
    if paired_input:
        log.info('Classifications were made from paired-end data.')
    else:
        log.info('Classifications were made from single-read data.')

    # Perform a naive check of the input file
    validate_input_file(args.original_classifications_file, verbose_input, args.minimum_hit_groups, paired_input)

    # Create a TaxonomyTree from the user provided names.dmp and nodes.dmp files
    taxonomy_tree = taxonomy.TaxonomyTree(names_filename=args.names, nodes_filename=args.nodes)

    # Filehandles-to-be
    o = None
    v = None

    # Open the classifications input file:
    with read_file(args.original_classifications_file) as f:
        log.info('Processing read classifications from "{file}".'.format(file=path.abspath(args.original_classifications_file)))

        # If user wants to save the read classifications to file, open file
        if args.output_classifications:
            if args.gz_output:
                if not args.output_classifications.endswith('.gz'):
                    args.output_classifications += '.gz'
            try:
                log.info('Saving reclassified reads in {}.'.format(args.output_classifications))
                o = write_file(args.output_classifications, args.gz_output)
            except (IOError, OSError, PermissionError) as e:
                log.error('Cannot open output classifications file "{}": {}'.format(args.output_classifications, e))
                sys.exit(1)

        # If user wants to save the verbose classification output to file, open file
        if args.output_verbose:
            if args.gz_output:
                if not args.output_verbose.endswith('.gz'):
                    args.output_verbose += '.gz'
            try:
                log.info('Saving verbose classification information in {}.'.format(args.output_verbose))
                v = write_file(args.output_verbose, args.gz_output)
            except (IOError, OSError, PermissionError) as e:
                log.error('Cannot open verbose output file "{}": {}'.format(args.output_verbose, e))
                sys.exit(1)

        # Run the main loop (reclassification)
        main_loop(f, tax_reads_dict, taxonomy_tree, args, report_frequency, taxa_lineages, paired_input, verbose_input, o, v)

    # Remember to close files
    if o:
        o.close()
    if v:
        v.close()


if __name__ == '__main__':
    stringmeup()
