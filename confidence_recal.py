#!/usr/bin/env python3

import argparse
import logging
import sys
import taxonomy
from os import path

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger(path.basename(__file__))

# TODO: supply more than one confidence threshold at the same time: $ python confidence_threshold.py 0.01 0.1 0.15 /path/to/file

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
    metavar='classifications_path',
    type=str,
    help='Path to the Kraken 2 output file containing the individual read classifications.')
parser.add_argument(
    '--names',
    metavar='taxonomy names'
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

taxonomy.TaxonomyTree()
