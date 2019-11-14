#!/usr/bin/env python3

import argparse
import logging
import sys
import taxonomy
from os import path

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger(path.basename(__file__))

# TODO: supply more than one confidence threshold at the same time: $ python confidence_threshold.py 0.01 0.1 0.15 /path/to/file
# TODO: make sure confidence_threshold is between 0 and 1
# TODO: be able to reclassify 'U' reads
# TODO: keep track of the number of reads at each taxon (after reclassification), so we can print a new report file more easily

# Arguments
parser = argparse.ArgumentParser(
    prog='Confidence score reclassification',
    usage='',
    description='Reclassification of reads (Kraken 2) by confidence score.')
parser.add_argument(
    'confidence_threshold',
    metavar='confidence',
    type=float,
    help='The confidence score threshold to be used in reclassification.')
parser.add_argument(
    'original_classifications_file',
    metavar='classifications_file',
    type=str,
    help='Path to the Kraken 2 output file containing the individual read classifications.')
parser.add_argument(
    '--names',
    metavar='taxonomy names'
)
parser.add_argument(
    '--nodes',
    metavar='taxonomy nodes'
)
parser.add_argument(
    '--output_report',
    type=str,
    action='store')
parser.add_argument(
    '--output_reads',
    type=str,
    action='store')
args = parser.parse_args()


def validate_input_file(putative_classifications_file):
    """
    Perform simple validation of the input file.
    """
    log.debug('Naive validation of input classifications file.')

    if not path.isfile(putative_classifications_file):
        log.error('Cannot find the classification file ({file}).'.format(
            file=putative_classifications_file))
        sys.exit()

    with open(putative_classifications_file, 'r') as f:
        line = f.readline()
        line_proc = line.strip()
        line_proc = line_proc.split('\t')

    # The following should be the case of a Kraken 2 output file
    num_cols = len(line_proc) == 5
    line_start = line_proc[0] in ['U', 'C']
    paired_data_1 = len(line_proc[3].split('|')) == 2
    paired_data_2 = len(line_proc[4].split('|:|')) == 2

    if num_cols and line_start and paired_data_1 and paired_data_2:
        log.debug('Validation OK.')
        return
    else:
        log.error('The classifications file is malformatted.')
        log.debug('First line of input: {line}'.format(line))
        log.debug('num_cols: {num_cols}'.format(num_cols))
        log.debug('line_start: {line_start}'.format(line_start))
        log.debug('paired_data_1: {paired_data_1}'.format(paired_data_1))
        log.debug('paired_data_2: {paired_data_2}'.format(paired_data_2))
        sys.exit()


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
    # Ambiguous hits are not processed (discarded).
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


def calc_confidence(current_tax_id, taxa_kmer_dict, taxonomy_tree):
    pass


def reclassify(classified_tax_id, taxa_kmer_dict, confidence_threshold, taxonomy_tree):
    """
    """
    current_node = classified_tax_id
    confidence_reached = False
    taxa_lineages = {}

    # The total number of assigned kmers (non-ambiguous):
    total_kmer_hits = sum(taxa_kmer_dict.values())

    # Only interested in tax_ids that are in the database. A '0' signifies that
    # the kmer could not be assigned to any tax_id.
    assigned_taxa_set = set([tax_id for tax_id in taxa_kmer_dict.keys() if tax_id != 0])
    print(assigned_taxa_set)
    while not confidence_reached:
        descendant_taxa = set()

        for tax_id in assigned_taxa_set:
            # TODO: Implement a 'is_descendant' fnc in taxonomy.py. Should return True as soon as the supplied tax_id is found in the upwards search. Should cut some time.
            lineage = taxa_lineages[tax_id] if taxa_lineages else taxonomy_tree.get_lineage([tax_id])[tax_id]
            if current_node in lineage:
                descendant_taxa.add(tax_id)

        # Sum the number of kmers that are assigned to any tax_id in the clade:
        num_hits_within_clade = sum([taxa_kmer_dict[tax_id] for tax_id in descendant_taxa])

        # The confidence value for the read pair classification at the current
        # taxonomic level:
        confidence = num_hits_within_clade / total_kmer_hits

        # If the confidence at this node is sufficient, we classify it to
        # the current node (TaxID).
        if confidence >= confidence_threshold:
            return (current_node, confidence)

        # Otherwise, set the current_node to the parent and keep going.
        else:
            current_node = taxonomy_tree.get_parent([current_node])[current_node]

            # If the next level is root, we shouldn't keep going.
            # The read pair is unclassified.
            if current_node == 1:
                return False

    # TODO: keep going up the taxonomy if the confidence score was not ok
    #print(confidence)
    #above_confidence =
    #while above_confidence:


if __name__ == '__main__':
    validate_input_file(args.original_classifications_file)
    taxonomy_tree = taxonomy.TaxonomyTree(names_filename=args.names, nodes_filename=args.nodes)

    with open(args.original_classifications_file, 'r') as f:
        log.info('Processing read classifications from "{file}".'.format(file=path.abspath(args.original_classifications_file)))
        i = 0
        for read_pair in f:
            read_pair_proc = read_pair.strip()
            read_pair_proc = read_pair_proc.split('\t')
            if read_pair_proc[0] == 'U':
                continue  # consider changing to read_pair.startswith('U') and put it first
            classified_tax_id = int(read_pair_proc[2])
            kmer_info_string = read_pair_proc[-1]
            taxa_kmer_dict = process_kmer_string(kmer_info_string)
            # TODO: taxa_kmer_dict can be empty if all kmer hits are ambiguous
            reclassify(classified_tax_id, taxa_kmer_dict, args.confidence_threshold, taxonomy_tree)
            i += 1
            if i % 1000000 == 0:
                print(i)
            # sys.exit()
