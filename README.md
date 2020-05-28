# StringMeUp

A post-processing tool to reclassify [Kraken 2] output based on the confidence score and/or minimum minimizer hit groups.

## Installation
StringMeUP is available to install through conda. Simply run the following command to install it:

`conda install -c conda-forge -c bioconda stringmeup`

## Usage

A good start is to run `stringmeup --help`.

## About the confidence score

The confidence score (CS) for a given read _R_ classified to a given node _J_ is calculated by dividing the number of k-mers that hit any node in the clade rooted at node _J_ (N) by the total number of k-mers that were queried against the database (M). Any k-mer with an ambiguous nucleotide is not queried against the database, and is thus not part of M.

CS = N / M

If the CS for a given read _R_ at a given node _J_ is equal to or larger than the specified cutoff, read _R_ is classified to node _J_. If not, the CS of read _R_ is calculated for the parent of node _J_. This is repeated until the CS >= CS cutoff or until we reach the root of the taxonomy. If the CS < CS cutoff at the root, the read is deemed unclassified.

## Reclassifying Kraken 2 output

To reclassify reads classified by Kraken 2 with a confidence cutoff of 0.1:

`stringmeup --names <names.dmp> --nodes <nodes.dmp> 0.1 <original_classifications.kraken2>`

Where:
* original_classifications.kraken2 is the output file from Kraken 2 that contain the read-by-read classifications.
* names.dmp and nodes.dmp are the same NCBI taxonomy files used for the building of the database that was used to produce the classifications in original_classifications.kraken2.

This command would output a Kraken 2 style report to stdout. Adding `--output_report <FILE>` would save the report in a file.

To save the read-by-read classifications, add `--output_classifications <FILE>` to the command.

To save a verbose version of the read-by-read classifications, add `--output_verbose <FILE>` to the command. The verbose version of the read-by-read classifications will contain the following columns:

| Column | Explanation |
|--------|-------------|
| READ_ID | The ID of the read |
| READ_LENGTH | The length of the read (same as Kraken 2 output) |
| MINIMIZER_HIT_GROUPS* | The number of minimizer hit groups found during Kraken 2 classification* |
| TAX_LVL_MOVES | How many levels in the taxonomy that the read moved during reclassification |
| ORIGINAL_TAXID | The taxID that the read was classified to originally |
| NEW_TAXID | The taxID that the read was reclassified to |
| ORIGINAL_CONFIDENCE | The original confidence score |
| NEW_CONFIDENCE | The confidence score at the taxID that the read was reclassified to |
| MAX_CONFIDENCE | The maximum confidence that the read can have |
| ORIGINAL_TAX_LVL | The taxonomic rank of the orignally classified taxID |
| NEW_TAX_LVL | The taxonomic rank of the reclassified taxID |
| ORIGINAL_NAME | The scientific name of the original taxID |
| NEW_NAME | The scientific name of the reclassified taxID |
| KMER_STRING | The k-mer string (same as Kraken 2 output) |

*: Is only present if the forked version of Kraken 2 was used for initial classification.

## Reclassifying with minimum hit groups

This option requires an input file that was produced with my [fork] of Kraken 2.

Add `--minimum_hit_groups <INT>` to the command. A read can only be considered classified if the number of minimizer hit groups is at or above the minimum_hit_groups setting.

[Kraken 2]: https://github.com/DerrickWood/kraken2
[fork]: https://github.com/danisven/kraken2
