#!/usr/bin/env python3

import argparse
import logging
import sys
import pickle
import datetime
from collections import namedtuple
from os import path

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger(path.basename(__file__))

# TODO: make it possible to use scientific names in the same way as tax_id

Node = namedtuple('Node', ['name', 'rank', 'parent', 'children'])

class TaxonomyTreeException(Exception):
    pass


class TaxonomyTree(object):
    """
    Creates a representation of the taxonomy in the files names.dmp and
    nodes.dmp of a kraken2 database.

    Inspired by https://github.com/frallain/NCBI_taxonomy_tree.
    """

    def __init__(self, nodes_filename=None, names_filename=None, pickled_taxonomy_filename=None):
        self.taxonomy = {}
        self.byranks = {}
        self.leaves = set()
        self.nodes_filename = nodes_filename
        self.names_filename = names_filename
        self.pickled_taxonomy_filename = pickled_taxonomy_filename
        #self.Node = namedtuple('Node', ['name', 'rank', 'parent', 'children'])

        if self.pickled_taxonomy_filename and (self.nodes_filename or self.names_filename):
            raise TaxonomyTreeException('You cannot combine nodes and names dump files with a pickled taxonomy file.')
        elif self.pickled_taxonomy_filename:
            self.read_pickled_taxonomy()
        elif self.nodes_filename and self.names_filename:
            self.construct_tree()
        else:
            pass

    def read_pickled_taxonomy(self, pickled_taxonomy_filename=None):
        """
        Reads a pickled taxonomy file and stores the taxonomy in self.taxonomy.
        """
        if not (self.pickled_taxonomy_filename or pickled_taxonomy_filename):
            raise TaxonomyTreeException('You tried to read a pickled taxonomy but didn\'t supply a filename.')
        pickle_file = pickled_taxonomy_filename if pickled_taxonomy_filename else self.pickled_taxonomy_filename

        try:
            log.info('Attempting to read pickled taxonomy from "{taxonomy_file}"...'.format(taxonomy_file=pickle_file))
            with open(pickle_file, 'rb') as f:
                taxonomy = pickle.load(f)
        except FileNotFoundError:
            log.exception('Could not find the pickled file "{taxonomy_file}".'.format(taxonomy_file=pickle_file))
            raise
        except:
            log.exception('There was an error when reading the pickled taxonomy file.')
            raise

        timestamp = taxonomy['timestamp']
        nodes_file = taxonomy['nodes_filename']
        names_file = taxonomy['names_filename']

        # TODO: perform some checks on the taxonomy to assert it is OK

        log.info('This taxonomy was created at timepoint {ts}.'.format(ts=timestamp.strftime("%Y-%m-%d %H:%M:%S")))
        log.info('This taxonomy was created from the nodes dump file "{nodes_filename}".'.format(nodes_filename=nodes_file))
        log.info('This taxonomy was created from the names dump file "{names_filename}".'.format(names_filename=names_file))
        self.taxonomy = taxonomy
        log.info('Done reading pickled taxonomy file.')

    def save_taxonomy(self, pickle_filename):
        """
        Stores the contents of self.taxonomy in a pickled file for later use.

        The pickled file can be read by an instance of TaxonomyTree instead of
        parsing names and nodes dump files.
        """
        log.info('Attempting to pickle the taxonomy to "{filename}".'.format(filename=pickle_filename))

        if not self.taxonomy:
            log.warning('NOTHING SAVED! There is no taxonomy to pickle. Returning.')
            return

        if path.isfile(pickle_filename):
            log.warning('NOTHING SAVED! A file with that name already exists ({filename}). Returning.'.format(filename=pickle_filename))
            return

        with open(pickle_filename, 'wb') as f:
            pickle.dump(self.taxonomy, f)

        log.info('Done pickling the taxonomy.')


    def construct_tree(self):
        """
        Reads a names.dmp and nodes.dmp file, and constructs a taxonomy tree
        representation:
            {tax_id#1: Node('name', 'rank', 'parent', 'children'),
             tax_id#2: Node('name', 'rank', 'parent', 'children'),
             ...,
             tax_id#N: ...}
        """

        if not (self.nodes_filename and self.names_filename):
            log.warning('You need to supply the nodes.dmp and names.dmp files.')
            return

        if self.taxonomy:
            log.warning('There was already a taxonomy tree, deleting the old one.')
            self.taxonomy = {}

        log.info("Constructing taxonomy tree...")
        taxid2name = {}

        try:
            log.info('Mapping taxonomic ID to scientific names from "{names_file}"...'.format(names_file=self.names_filename))
            # TODO: check so that names.dmp conforms to expected format
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
            log.info('Reading taxonomy from "{nodes_file}"...'.format(nodes_file=self.nodes_filename))
            # TODO: check so that nodes.dmp conforms to expected format
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
                        node = Node(name=tax_name, rank=tax_rank, parent=tax_parent, children=[])
                        self.taxonomy[tax_id] = node
                        self.leaves.add(tax_id)

                    if tax_parent in self.taxonomy:
                        self.taxonomy[tax_parent].children.append(tax_id)
                        if tax_parent in self.leaves:
                            self.leaves.remove(tax_parent)
                    else:
                        parent_node = Node(name=taxid2name[tax_parent], rank=None, parent=None, children=[tax_id])
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
        self.taxonomy['timestamp'] = datetime.datetime.now()
        self.taxonomy['nodes_filename'] = path.abspath(self.nodes_filename)
        self.taxonomy['names_filename'] = path.abspath(self.names_filename)
        log.info("Taxonomy tree built.")

    def _get_property(self, tax_id, property):
        """
        Internal function to fetch the value of a single property of a namedtuple in the taxonomy dictionary.
        Raises an exception if tax_id does not exist in the taxonomy tree.
        Raises an exception if the taxonomy tree isn't built yet.
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
            log.exception('You have not built the taxonomy tree yet.')
            raise TaxonomyTreeException('You have not built the taxonomy tree yet.')

        return property_value

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

    def get_node(self, tax_id_list):
        """
        Returns the node instances of the supplied tax_ids.
        """
        #TODO: Use this fnc in other fncs when getting nodes from self.taxonomy
        self._verify_list(tax_id_list)
        node_dict = {}

        if self.taxonomy:
            for tax_id in tax_id_list:
                try:
                    node = self.taxonomy[tax_id]
                except KeyError:
                    log.exception('Could not find tax_id={tax_id} in the taxonomy tree.'.format(tax_id=tax_id))
                    raise
                node_dict[tax_id] = node
        else:
            log.exception('You have not built the taxonomy tree yet.')
            raise TaxonomyTreeException('You have not built the taxonomy tree yet.')

        return node_dict

    def get_lineage(self, tax_id_list):
        """
        For each tax_id, returns the input tax_id and the tax_ids of its
        ancestors.
        """
        self._verify_list(tax_id_list)
        lineage_dict = {}

        for tax_id in tax_id_list:
            lineage = [tax_id]
            node = self.get_node([tax_id])[tax_id]

            while node.parent:
                lineage.append(node.parent)
                node = self.get_node([node.parent])[node.parent]

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
            node = self.get_node([tax_id])[tax_id]
            children_pool = set(node.children)
            clade = set([tax_id])
            clade.update(children_pool)

            while children_pool:
                try:
                    clade_taxon = children_pool.pop()
                except KeyError:
                    break
                else:
                    new_children = self.get_node([clade_taxon])[clade_taxon].children
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

    def get_clade_rank_taxids(self, tax_ids, rank):
        """
        For each clade rooted at the input tax_ids, return all tax_ids that
        represent taxa at the supplied rank. For example:
        # get_clade_rank_taxids([1], 'phylum') -- returns all phyla in the whole tree
        # get_clade_rank_taxids([2, 9443], 'genus') -- returns all genera in the clades rooted at 'Bacteria' and 'Primates'
        """
        self._verify_list(tax_ids)
        clade_dict = {}
        for tax_id in tax_ids:
            clade = self.get_clade([tax_id])[tax_id]
            clade_rank_taxids = None

            try:
                clade_rank_taxids = clade.intersection(self.byranks[rank])
            except KeyError:
                log.warning('No such rank: {rank}.'.format(rank=rank))

            clade_dict[tax_id] = clade_rank_taxids

        return clade_dict

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