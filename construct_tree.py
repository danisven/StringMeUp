#!/usr/bin/env python3

import argparse
import logging
import sys
from collections import namedtuple
from os import path

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger(path.basename(__file__))


class TaxonomyTree(object):
    """
    Creates a representation of the taxonomy in the files names.dmp and
    nodes.dmp of a kraken2 database.
    """

    def __init__(self, nodes_filename=None, names_filename=None):
        self.taxonomy = {}
        self.nodes_filename = nodes_filename
        self.names_filename = names_filename
        self.Node = namedtuple('Node', ['name', 'rank', 'parent', 'children'])

        if self.names_filename and self.nodes_filename:
            log.info("Constructing taxonomy tree...")
            self._construct_taxonomy_tree()


    def _construct_taxonomy_tree(self):
        log.info("Mapping taxonomic ID to scientific names from {names_file}...".format(names_file=self.names_filename))
        taxid2name = {}
        with open(self.names_filename, 'r') as f:
            for name_line in f:
                name_info = name_line.split('|')
                name_type = name_info[3].strip()
                if name_type == 'scientific name':
                    tax_id = int(name_info[0].strip())
                    tax_name = name_info[1].strip()
                    taxid2name[tax_id] = tax_name

        log.info("Reading taxonomy from {nodes_file}...".format(nodes_file=self.nodes_filename))
        with open(self.nodes_filename, 'r') as f:
            for tax_line in f:
                tax_info = tax_line.split('|')[0:3]
                tax_id = int(tax_info[0].strip())
                tax_parent = int(tax_info[1].strip())
                tax_rank = tax_info[2].strip()
                tax_name = taxid2name[tax_id]

                if tax_id in self.taxonomy:
                    # We already inserted the current tax_id as a parent of another
                    self.taxonomy[tax_id] = self.taxonomy[tax_id]._replace(rank=tax_rank, parent=tax_parent)
                else:
                    node = self.Node(name=tax_name, rank=tax_rank, parent=tax_parent, children=[])
                    self.taxonomy[tax_id] = node

                if tax_parent in self.taxonomy:
                    self.taxonomy[tax_parent].children.append(tax_id)
                else:
                    parent_node = self.Node(name=taxid2name[tax_parent], rank=None, parent=None, children=[tax_id])
                    self.taxonomy[tax_parent] = parent_node


    def _get_property(self, tax_id, property):
        """
        Internal function to fetch the value of a single property of a namedtuple in the taxonomy dictionary.
        """

        if self.taxonomy:
            try:
                property_value = getattr(self.taxonomy[tax_id], property)
            except KeyError as error:
                log.exception('Could not find tax_id={tax_id} in the taxonomy tree.'.format(tax_id=tax_id), err.message)
                raise
            except AttributeError as error:
                log.exception('There is no such field ("{field}") in the namedtuple.'.format(field=property), err.message)
                raise
            else:
                return property_value
        else:
            #TODO: handle
            raise


    def get_children(self, tax_id):
        """
        Returns the direct descending children of a tax_id (one level).
        """

        return self._get_property(tax_id, 'children')


    def get_parent(self, tax_id):
        """
        Returns the parent of a tax_id.
        """

        return self._get_property(tax_id, 'parent')
        # TODO: what happens at the root?


    def get_rank(self, tax_id):
        """
        Returns the rank of a tax_id.
        """

        return self._get_property(tax_id, 'rank')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--nodes')
    parser.add_argument('--names')
    args = parser.parse_args()

    taxonomy_tree = TaxonomyTree(args.nodes, args.names)
