#!/usr/bin/env python3

__version__ = "0.1.1"

import argparse
import stringmeup.taxonomy
import operator
import logging
import gzip
import sys
from dataclasses import dataclass
from os import path

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


def validate_input_file(putative_classifications_file, verbose_input, minimum_hit_groups):
    """
    Perform simple validation of the input file.
    """
    log.debug('Naive validation of input classifications file.')

    if not path.isfile(putative_classifications_file):
        log.error('Cannot find the classification file ({file}).'.format(
            file=putative_classifications_file))
        sys.exit()

    with read_file(putative_classifications_file) as f:
        line = f.readline()
        line_proc = line.strip()
        line_proc = line_proc.split('\t')

        # The following should be the case of a Kraken 2 output file
        # First, check so the number of columns in the input file conforms to the expected number (especially when using/not using minimum_hit_groups):
        # TODO: This can be done much better.
        if not verbose_input:
            num_cols = len(line_proc) == 5  # original type of kraken2 output file
            if minimum_hit_groups:
                log.error('You specified --minimum_hit_groups {}, but didn\'t supply an input file that contain a minimizer hit groups column.')
                sys.exit()
        else:
            num_cols = len(line_proc) == 6  # 6 columns if the output was produced with the verbose version of kraken2 that outputs minimizer_groups
            if not minimum_hit_groups:
                log.error('The input file contains too many columns. Use --minimum_hit_groups if you want to use an input file that contain a minimizer hit groups column.')
                sys.exit()

        line_start = line_proc[0] in ['U', 'C']
        paired_data_1 = len(line_proc[3].split('|')) == 2
        paired_data_2 = len(line_proc[-1].split('|:|')) == 2  # Should be enough to change this line if we want to accomodate reclassification of single reads

        if num_cols and line_start and paired_data_1 and paired_data_2:
            log.debug('Validation OK.')
            return
        else:
            log.error('The classifications file is malformatted.')
            log.debug('First line of input: {}'.format(line))
            log.debug('num_cols: {}'.format(num_cols))
            log.debug('line_start: {}'.format(line_start))
            log.debug('paired_data_1: {}'.format(paired_data_1))
            log.debug('paired_data_2: {}'.format(paired_data_2))
            sys.exit()


def is_verbose_input(classifications_file):
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


def process_kmer_string(kmer_info_string):
    """
    Process a kmer info string (last column of a Kraken 2 output file), so that
    we get a dictionary mapping of tax_ids to total sum of kmer hits.
    Returns:
    {tax_id_#1: X kmer hits,
     tax_id_#2: Y kmer hits,
     ...
     tax_id_#N: Z kmer hits}
    """
    kmer_info_string = kmer_info_string.split()
    kmer_info_string.remove('|:|')

    # Messy list comprehension. Converts all "taxa":"num_kmer" string pairs
    # into integer tuples like (taxa, num_kmers), and saves them in a list.
    # Ambiguous kmers are not processed (discarded).
    kmer_classifications = [
        (int(x[0]), int(x[1])) for x in (
            kmer_info.split(':') for kmer_info in kmer_info_string)
        if x[0] != 'A']

    # Further processes the (taxa, num_kmers) tuples into a dict where each
    # tax_id stores the total sum of kmer hits to that tax_id.
    taxa_kmer_dict = {}
    for kmer_info in kmer_classifications:
        if kmer_info[0] not in taxa_kmer_dict:
            taxa_kmer_dict[kmer_info[0]] = kmer_info[1]
        else:
            taxa_kmer_dict[kmer_info[0]] += kmer_info[1]

    return taxa_kmer_dict


def reclassify_read(read, confidence_threshold, taxonomy_tree, verbose_input, minimum_hit_groups, taxa_lineages):
    """
    Sums the number of kmers that hit in the clade rooted at "current_node",
    and divides it with the total number of kmers queried against the database:
    confidence = clade_kmer_hits / total_kmer_hits

    If the confidence at a specific node is < confidence_threshold, we go one
    step up the taxonomy (to the parent node) and recalculates the confidence.
    This is repeated until confidence >= confidence_threshold.

    In this function it's envisionable to include other parameters for the
    classification... Right now I'm only considering the confidence score.
    """
    # Process the kmer string into a dict of {tax_id: #kmers} key, value pairs
    taxa_kmer_dict = process_kmer_string(read.kmer_string)

    # Make the current node the same as the original classification
    read.current_node = read.original_taxid

    # The total number of kmers that were interrogated against the
    # database (non-ambiguous):
    total_kmer_hits = sum(taxa_kmer_dict.values())

    # Only interested in tax_ids that are in the database. A '0' signifies that
    # the kmer could not be assigned to any tax_id (missing from database).
    assigned_taxa_set = set(
        [tax_id for tax_id in taxa_kmer_dict.keys() if tax_id != 0])

    # Make a quick check to see if it is even possible to obtain the confidence
    # needed to make a classification. If it isn't we don't have to go through
    # the hassle of calculating the confidence at all parent nodes. Potentially
    # saving us a lot of time.
    doomed_to_fail = False
    total_hits = sum(taxa_kmer_dict[tax_id] for tax_id in assigned_taxa_set)
    max_confidence = total_hits / total_kmer_hits
    read.max_confidence = max_confidence

    # The read can't achieve a confidence high enough, so we mark it
    if max_confidence < confidence_threshold:
        doomed_to_fail = True

    # Filter minimizer_hit_groups
    if verbose_input:
        if read.minimizer_hit_groups < minimum_hit_groups:
            doomed_to_fail = True

    # The nr of kmers that hit within the clade rooted at the current node:
    num_hits_within_clade = 0

    while not read.classified:
        taxa_in_clade = set()

        # For each tax_id that kmers in the read were assigned to:
        for tax_id in assigned_taxa_set:

            # Get the lineage (all ancestors including itself) for the tax_id:
            if tax_id in taxa_lineages:
                lineage = taxa_lineages[tax_id]
            else:
                lineage = taxonomy_tree.get_lineage([tax_id])[tax_id]

                # Save lineage so we don't have to get it from taxonomy_tree
                # more than once. Also make it into a set, which is faster to
                # query (the order of tax_ids in the lineage is not important
                # here).
                taxa_lineages[tax_id] = set(lineage)

            # If the currently classified (read.current_node) tax_id is in the
            # lineage (parents) of tax_id, then tax_id must be in the clade
            # rooted at read.current_node - i.e. tax_id is a descendant of
            # current_node.
            if read.current_node in lineage:

                # There is no need to get the lineage of tax_id in future
                # iterations since it will always be in the clade rooted at
                # read.current_node (we only ever go up in the taxonomy).
                # Remember which tax_ids we have counted, so we can remove them
                # from the set, outside of the loop:
                taxa_in_clade.add(tax_id)

                # Instead, we just add the kmers that hit tax_id to the total
                # hits at the clade:
                num_hits_within_clade += taxa_kmer_dict[tax_id]

        # Remove the already counted tax_ids:
        if taxa_in_clade:
            assigned_taxa_set -= taxa_in_clade

        # The confidence value for the read pair classification at the current
        # taxonomic level:
        read.recalculated_conf = num_hits_within_clade / total_kmer_hits

        # Set the original confidence score
        if not read.original_conf:
            read.original_conf = read.recalculated_conf

            # If we can't achieve the confidence score cutoff, now is the time
            # to exit the loop and return the read (since we have calculated
            # the original confidence).
            if doomed_to_fail:
                read.recalculated_conf = max_confidence
                read.reclassified_taxid = 0
                return read, taxa_lineages

        # If the confidence at this node is sufficient, we classify it to
        # the current node (TaxID).
        if read.recalculated_conf >= confidence_threshold:
            read.classified = True
            read.reclassified_taxid = read.current_node
            return read, taxa_lineages

        # If the current node is the root, we stop (can't go higher up in the
        # taxonomy).
        elif read.current_node == 1:
            read.reclassified_taxid = 0
            return read, taxa_lineages

        # Otherwise, set the current_node to the parent and keep going.
        else:
            read.current_node = taxonomy_tree.get_parent(
                [read.current_node])[read.current_node]


def read_kraken_output():
    """
    This should work to read classifications from a kraken 2 output file. It's
    not complete, but the backbone is there. The point is, we should not have
    to run reclassification in order to produce a report file - we shuold be
    able to just read an output file and work with the classifications as they
    are. Could also just modify the main loop.
    """
    tax_dict = {'hits_at_node': {}, 'hits_at_clade': {}}
    i = 0
    with open('Ki-1974-23-291.kraken2', 'r') as f:
        for line in f:
            if line.startswith('C'):
                l = line.strip().split('\t')
                tax_id = int(l[2])
                if tax_id not in tax_dict['hits_at_node']:
                    tax_dict['hits_at_node'][tax_id] = 1
                else:
                    tax_dict['hits_at_node'][tax_id] += 1
            i += 1
            if i % report_frequency == 0:
                print('Processed {} lines'.format(i))


def get_kraken2_report_content(tax_reads, taxonomy_tree, total_reads):
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


def format_kraken2_report_row(report_node):
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


def make_kraken2_report(tax_reads, taxonomy_tree, total_reads, output_report):
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


def get_verbose_output(read, taxonomy_tree):
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


def create_read(kraken2_read, verbose_input=False):
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


def main_loop(f_handle, tax_reads_dict, taxonomy_tree, args, report_frequency, taxa_lineages, verbose_input=False, o_handle=None, v_handle=None):
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

            # Reclassify the read pair based on confidence
            read, taxa_lineages = reclassify_read(
                read,
                args.confidence_threshold,
                taxonomy_tree,
                verbose_input,
                args.minimum_hit_groups,
                taxa_lineages)

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
        if i % report_frequency == 0:
            log.info('Processed {} reads...'.format(i))

    log.info('Done processing reads. They were {} in total.'.format(i))

    # Output a report file
    make_kraken2_report(tax_reads_dict, taxonomy_tree, i, args.output_report)  # i is used to calculate the ratio of classified reads (col 1 in output file).


def read_file(filename):
    """
    Wrapper to read either gzipped or ordinary text file input.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


def write_file(filename, gz_output):
    """
    Wrapper to write either gzipped or ordinary text file output.
    """
    if gz_output:
        return gzip.open(filename, 'wt')
    else:
        return open(filename, 'w')


def get_arguments():
    """
    Wrapper fnc to get the command line arguments. Inserting this piece of code
    into its own fnc for conda compatibility.
    """

    parser = argparse.ArgumentParser(
        prog='Confidence score reclassification',
        usage='confidence_recal [--names <FILE> --nodes <FILE> | --taxonomy_file <FILE>] [--output_report <FILE>] [--output_classifications <FILE>] [--output_verbose <FILE>] [--keep_unclassified] [--minimum_hit_groups INT] [--gz_output] [--help] confidence classifications',
        description='Reclassification of reads (Kraken 2) by confidence score.')
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
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--names',
        metavar='FILE',
        help='taxonomy names dump file (names.dmp)')
    group.add_argument(
        '--taxonomy_file',
        metavar='FILE',
        help='Path to a pickle file containing a taxonomy tree created through the TaxonomyTree.save_taxonomy function (taxonomy.py).')
    parser.add_argument(
        '--nodes',
        metavar='FILE',
        help='taxonomy nodes dump file (nodes.dmp)')
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


def stringmeup():

    # Get the CL arguments
    args = get_arguments()

    if args.names and args.nodes is None:
        parser.error("--names requires --nodes to be set as well.")

    # Some initial setup
    taxa_lineages = {}
    report_frequency = 10000000  # Will output progress every nth read
    tax_reads_dict = {'hits_at_node': {}, 'hits_at_clade': {}}

    # Was the input generated with https://github.com/danisven/kraken2 ?
    verbose_input = is_verbose_input(args.original_classifications_file)
    if verbose_input:
        log.info("The input file appears to contain minimizer_hit_groups.")

    # Perform a naive check of the input file
    validate_input_file(args.original_classifications_file, verbose_input, args.minimum_hit_groups)

    # If user provided names.dmp and nodes.dmp, create taxonomy tree from that,
    # otherwise, create i from a pickled taxonomy file
    if args.names:
        taxonomy_tree = taxonomy.TaxonomyTree(names_filename=args.names, nodes_filename=args.nodes)
    else:
        taxonomy_tree = taxonomy.TaxonomyTree(pickled_taxonomy_filename=args.taxonomy_file)

    # Filehandles-to-be
    o = None
    v = None

    # Open the classifications input file:
    with read_file(args.original_classifications_file) as f:
        log.info('Processing read classifications from "{file}".'.format(file=path.abspath(args.original_classifications_file)))

        # TODO: make sure output files are writable
        # If user wants to save the read classifications to file, open file
        if args.output_classifications:
            if args.gz_output:
                if not args.output_classifications.endswith('.gz'):
                    args.output_classifications += '.gz'
            log.info('Saving reclassified reads in {}.'.format(args.output_classifications))
            o = write_file(args.output_classifications, args.gz_output)

        # If user wants to save the verbose classification output to file, open file
        if args.output_verbose:
            if args.gz_output:
                if not args.output_verbose.endswith('.gz'):
                    args.output_verbose += '.gz'
            log.info('Saving verbose classification information in {}.'.format(args.output_verbose))
            v = write_file(args.output_verbose, args.gz_output)

        # Run the main loop (reclassification)
        main_loop(f, tax_reads_dict, taxonomy_tree, args, report_frequency, taxa_lineages, verbose_input, o, v)

    # Remember to close files
    if o:
        o.close()
    if v:
        v.close()


if __name__ == '__main__':
    stringmeup()
