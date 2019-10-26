#!/usr/bin/env python3

import argparse
import logging
import sys
from collections import namedtuple
from os import path

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger(path.basename(__file__))

# TODO: make it possible to use scientific names in the same way as tax_id

class TaxonomyTreeException(Exception):
    pass


class TaxonomyTree(object):
    """
    Creates a representation of the taxonomy in the files names.dmp and
    nodes.dmp of a kraken2 database.

    Inspired by https://github.com/frallain/NCBI_taxonomy_tree.
    """

    def __init__(self, nodes_filename=None, names_filename=None):
        self.taxonomy = {}
        self.byranks = {}
        self.leaves = set()
        self.nodes_filename = nodes_filename
        self.names_filename = names_filename
        self.Node = namedtuple('Node', ['name', 'rank', 'parent', 'children'])

        if self.names_filename and self.nodes_filename:
            self.construct_tree()

    def construct_tree(self):
        """
        Reads a names.dmp and nodes.dmp file, and constructs a taxonomy tree
        representation:
            {tax_id#1: Node('name', 'rank', 'parent', 'children'),
             tax_id#2: Node('name', 'rank', 'parent', 'children'),
             ...,
             tax_id#N: ...}
        """
        if self.taxonomy:
            log.info('There was already a taxonomy tree, deleting the old one.')
            self.taxonomy = {}

        log.info("Constructing taxonomy tree...")
        taxid2name = {}

        try:
            log.info("Mapping taxonomic ID to scientific names from {names_file}...".format(names_file=self.names_filename))
            with open(self.names_filename, 'r') as f:
                for name_line in f:
                    name_info = name_line.split('|')
                    name_type = name_info[3].strip()
                    if name_type == 'scientific name':
                        tax_id = int(name_info[0].strip())
                        tax_name = name_info[1].strip()
                        taxid2name[tax_id] = tax_name
        except FileNotFoundError:
            log.exception('Could not find the file "{names_file}".'.format(names_file=self.names_filename))
            raise

        try:
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
                        self.leaves.add(tax_id)

                    if tax_parent in self.taxonomy:
                        self.taxonomy[tax_parent].children.append(tax_id)
                        if tax_parent in self.leaves:
                            self.leaves.remove(tax_parent)
                    else:
                        parent_node = self.Node(name=taxid2name[tax_parent], rank=None, parent=None, children=[tax_id])
                        self.taxonomy[tax_parent] = parent_node

                    # Save the tax_id to it's corresponding rank set
                    if tax_rank in self.byranks:
                        self.byranks[tax_rank].add(tax_id)
                    else:
                        self.byranks[tax_rank] = set([tax_id])

        except FileNotFoundError:
            log.exception('Could not find the nodes file "{nodes_file}".'.format(nodes_file=self.nodes_filename))
            raise

        # Adjust the root
        root_children = self.taxonomy[1].children
        root_children.remove(1)
        self.taxonomy[1] = self.taxonomy[1]._replace(parent=None, children=root_children)
        log.info("Taxonomy tree built.")

    def _get_property(self, tax_id, property):
        """
        Internal function to fetch the value of a single property of a namedtuple in the taxonomy dictionary.
        """
        if self.taxonomy:
            try:
                property_value = getattr(self.taxonomy[tax_id], property)
            except KeyError:
                log.exception('Could not find tax_id={tax_id} in the taxonomy tree.'.format(tax_id=tax_id))
                raise
            except AttributeError:
                log.exception('There is no such field ("{field}") in the namedtuple.'.format(field=property))
                raise
            else:
                return property_value
        else:
            log.exception('You have not built the taxonomy tree yet.')
            raise TaxonomyTreeException('You have not built the taxonomy tree yet.')

    def _verify_list(self, putative_list):
        """
        Internal helper function to check that input lists are indeed lists.
        """
        try:
            assert isinstance(putative_list, list)
        except AssertionError:
            log.exception('Input must be a list. You input "{input}", of type {input_type}'.format(
                input=putative_list, input_type=type(putative_list)))
            raise

    def get_name(self, tax_id_list):
        """
        Returns the names of the tax_ids in the input list.
        """
        self._verify_list(tax_id_list)
        name_dict = {}
        for tax_id in tax_id_list:
            name_dict[tax_id] = self._get_property(tax_id, 'name')
        return name_dict

    def get_children(self, tax_id_list):
        """
        Returns the direct descending children of each tax_id.
        """
        self._verify_list(tax_id_list)
        children_dict = {}
        for tax_id in tax_id_list:
            children_dict[tax_id] = self._get_property(tax_id, 'children')
        return children_dict

    def get_parent(self, tax_id_list):
        """
        Returns the parent of each tax_id.
        """
        self._verify_list(tax_id_list)
        parent_dict = {}
        for tax_id in tax_id_list:
            parent_dict[tax_id] = self._get_property(tax_id, 'parent')
        return parent_dict

    def get_rank(self, tax_id_list):
        """
        Returns the rank of each tax_id.
        """
        self._verify_list(tax_id_list)
        rank_dict = {}
        for tax_id in tax_id_list:
            rank_dict[tax_id] = self._get_property(tax_id, 'rank')
        return rank_dict

    def get_lineage(self, tax_id_list):
        """
        For each tax_id, returns the input tax_id and the tax_ids of its
        ancestors.
        """
        self._verify_list(tax_id_list)
        lineage_dict = {}

        for tax_id in tax_id_list:
            lineage = [tax_id]
            node = self.taxonomy[tax_id]

            while node.parent:
                lineage.append(node.parent)
                node = self.taxonomy[node.parent]

            lineage.reverse()
            lineage_dict[tax_id] = lineage

        return lineage_dict

    def get_clade(self, tax_id_list):
        """
        For each tax_id, returns all of the tax_ids of the clade rooted at the
        tax_id.

        returns: {tax_id#1: set(all tax_ids in node),
                  tax_id#2: set(all tax_ids in node)}
        """

        self._verify_list(tax_id_list)
        clade_dict = {}

        for tax_id in tax_id_list:
            node = self.taxonomy[tax_id]
            children_pool = set(node.children)
            clade = set([tax_id])
            clade.update(children_pool)

            while children_pool:
                try:
                    clade_taxon = children_pool.pop()
                except KeyError:
                    break
                else:
                    new_children = self.taxonomy[clade_taxon].children
                    clade.update(new_children)
                    children_pool.update(new_children)

            clade_dict[tax_id] = clade

        return clade_dict

    def get_leaves(self, tax_ids=[1]):
        """
        Returns a {tax_id: set(leaf_taxids)} mapping of leaf node tax_ids for
        the clades rooted at the tax_ids.
        """

        self._verify_list(tax_ids)
        clade_dict = {}
        if tax_ids == [1]:
            clade_dict[1] = self.leaves

        for tax_id in tax_ids:
            clade = self.get_clade([tax_id])[tax_id]
            clade_leaves = self.leaves.intersection(clade)
            clade_dict[tax_id] = clade_leaves

        return clade_dict

    def get_clade_species(self, tax_ids):
        """
        For each clade rooted at the input tax_ids, return the tax_ids of the
        species of that clade.
        """
        # TODO: must verify that the supplied tax_ids are at or above the rank of species
        pass

    def get_clade_rank_taxids(self, tax_ids, rank):
        """
        For each clade rooted at the input tax_ids, return all tax_ids that
        represent taxa at the supplied rank. For example:
        # get_clade_rank_taxids([1], 'phylum') -- returns all phyla in the whole tree
        # get_clade_rank_taxids([2, 9443], 'genus') -- returns all genera in the clades rooted at 'Bacteria' and 'Primates'
        """
        # TODO: check that the supplied tax_ids are at or above the supplied rank. Should log a warning if not.
        return self.byranks[rank].intersection(self.get_clade(tax_ids)[tax_ids[0]])

    def set_taxonomy_files(self, nodes_filename, names_filename):
        log.info('Setting the nodes filename to {}'.format(nodes_filename))
        self.nodes_filename = nodes_filename
        log.info('Setting the names filename to {}'.format(names_filename))
        self.names_filename = names_filename

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--nodes')
    parser.add_argument('--names')
    args = parser.parse_args()

    taxonomy_tree = TaxonomyTree(args.nodes, args.names)
