"""
These methods are used for grouping CoMetGeNe trails by reactions.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import sys

import numpy as np

from NeighborhoodResolver import NeighborhoodResolver
from ..definitions import PICKLE_GEN
from trail_grouping_utils import get_species_object
from trail.utils import split_entities, pickle
from phylogeny import get_ordered_organisms


def add_to_species_dict(species, species_dict):
    """Builds genome information for the specified species, and stores it in a
    pickle file.

    :param species: species for which to build genome information
    :param species_dict: dict with genome information for other species
    """
    species_dict[species] = get_species_object(species)
    pickle(PICKLE_GEN, species_dict)


def get_reactions_in_sets(species_reaction_sets):
    """Determines and returns the list of reactions in the reaction sets for
    the query species.

    :param species_reaction_sets: list of species-specific CoMetGeNe reaction
        sets
    :return: list of reactions in species_reaction_sets
    """
    reactions_in_sets = list()

    for reaction_set in species_reaction_sets:
        for reaction in reaction_set.reactions:
            reactions_in_sets.append(reaction)

    return reactions_in_sets


def get_symbols(init_ds, species, kgml_dir):
    """Determines symbols for the reaction table of species.

    :param init_ds: list of necessary data structures storing R number
        associations, general and species-specific CoMetGeNe reaction sets, and
        genome information for all species
    :param species: KEGG organism code for the reference species
    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :return: numpy character array with symbols
    """
    assert len(init_ds) == 4
    r_numbers = init_ds[0]
    reaction_sets = init_ds[1]
    species_reaction_sets = init_ds[2]
    species_dict = init_ds[3]

    reactions_in_sets = get_reactions_in_sets(species_reaction_sets)
    organisms = get_ordered_organisms(species, kgml_dir)
    symbols = np.chararray((len(reactions_in_sets), len(organisms)))

    for srs_idx in range(len(species_reaction_sets)):
        reaction_set = species_reaction_sets[srs_idx]
        rs = reaction_sets.find(reaction_set.reactions)

        ris_idx = 0  # index of the 1st reaction in rs in reactions_in_sets
        for i in range(srs_idx):
            r_set = species_reaction_sets[i]
            ris_idx += len(r_set.reactions)

        for org_idx in range(len(organisms)):
            organism = organisms[org_idx]
            if organism not in species_dict:
                add_to_species_dict(organism, species_dict)
            other_species = species_dict[organism]
            n_resolver = NeighborhoodResolver(other_species, rs, r_numbers)
            for i in range(len(rs.reactions)):
                reaction = reactions_in_sets[ris_idx + i]
                symbols[ris_idx + i, org_idx] = n_resolver.state(reaction)

    return symbols


def get_genes_and_pathways(reactions, r_numbers, species):
    """Returns a CSV-formatted string with the list of genes and pathways where
    the reaction(s) of 'species' appear.

    :param reactions: list of reactions for species
    :param r_numbers: RNumbers object
    :param species: KEGG organism code
    :return: CSV-formatted string with genes and pathways where reactions of
        species are present
    """
    gene_set = set()
    pathway_set = set()

    for reaction in reactions:
        organism = r_numbers.find(reaction).find(species)
        assert organism is not None
        for gene in organism.genes:
            gene_set.add(gene.replace(species + ':', ''))
        for pathway in organism.pathways:
            pathway_set.add(pathway)

    gene_col = ' '.join(sorted(gene_set))
    pathway_col = ' '.join(sorted(pathway_set))

    return gene_col.rstrip() + ';' + pathway_col.rstrip() + ';'


def find_reaction_set_index(i, species_reaction_sets):
    """Returns the index of the reaction set in species_reaction_sets where the
    i-th reaction is found.

    :param i: index of a reaction in a reaction set from species_reaction_sets
    :param species_reaction_sets: list of species-specific CoMetGeNe reaction
        sets
    :return: index in species_reaction_sets for the reaction
    """
    nb_reactions = 0
    srs_index = 0

    while nb_reactions < i:
        reaction_set = species_reaction_sets[srs_index]
        nb_reactions += len(reaction_set.reactions)
        srs_index += 1

    return srs_index - 1


def find_last_reaction_index(species_reaction_sets, srs_idx):
    """Returns the index in species_reaction_sets where the reaction set at
    srs_idx ends.

    :param species_reaction_sets: list of species-specific CoMetGeNe reaction
        sets
    :param srs_idx: index of a reaction set in species_reaction_sets
    :return: index in terms of reactions in species_reaction_sets where the
        reaction set at index srs_idx ends
    """
    last_idx = 0

    for i in range(srs_idx + 1):
        reaction_set = species_reaction_sets[i]
        last_idx += len(reaction_set.reactions)

    return last_idx - 1


def print_header(species, organisms, dev_out):
    """Prints table header for the reference species to the specified output
    device.

    :param species: KEGG organism code for the reference species
    :param organisms: list of other species in the data set
    :param dev_out: output device (file handle or stdout)
    """
    header = 'reaction;' + species + '_gene;pathway;'
    for i in range(len(organisms)):
        header += organisms[i]
        if i < len(organisms) - 1:
            header += ';'
    dev_out.write(header + '\n')


def print_symbols(init_ds, symbols, species, kgml_dir, dev_out):
    """Prints the table 'symbols' representing trail grouping by reactions for
    the reference species to the output device dev_out.

    :param init_ds: list of necessary data structures storing R number
        associations, general and species-specific CoMetGeNe reaction sets, and
        genome information for all species
    :param symbols: numpy character array with symbols to print to dev_out
    :param species: KEGG organism code for the reference species
    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :param dev_out: output device (file handle or stdout)
    """
    organisms = get_ordered_organisms(species, kgml_dir)
    print_header(species, organisms, dev_out)

    assert len(init_ds) == 4
    r_numbers = init_ds[0]
    species_reaction_sets = init_ds[2]
    reactions_in_sets = get_reactions_in_sets(species_reaction_sets)

    for i in range(symbols.shape[0]):
        line = reactions_in_sets[i] + ';'
        reactions = split_entities([reactions_in_sets[i]])
        line += get_genes_and_pathways(reactions, r_numbers, species)

        for j in range(symbols.shape[1]):
            line += symbols[i, j]
            if j < symbols.shape[1] - 1:
                line += ';'
        dev_out.write(line + '\n')

        srs_idx = find_reaction_set_index(i, species_reaction_sets)
        end_rs_index = find_last_reaction_index(species_reaction_sets, srs_idx)
        if i == end_rs_index:
            dev_out.write('***\n')

    if dev_out != sys.stdout:
        dev_out.close()
