"""
These utilities are needed for grouping CoMetGeNe trails by genes and/or by
reactions.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import os
import sys

from ReactionSet import ReactionSets
from Species import Chromosome, Species
from RNumber import RNumbers
from TrailGroupingError import TrailGroupingError
from ..definitions import PICKLE_GENOME, PICKLE_RN, PICKLE_RS, PICKLE_GEN
from trail.utils import pickle, unpickle
from parsers.kgml_parser import parse_kgml


def initialization(results, species, kgml_dir, by_r):
    """Initializes and returns data structures required for trail grouping.

    :param results: parsed CoMetGeNe results
    :param species: reference species (KEGG organism code)
    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :param by_r: True if grouping by reactions, False if grouping by genes
    :return: list of data structures storing R number associations, CoMetGeNe
        reaction sets, and genome information for all species (as well as
        species-specific CoMetGeNe reaction sets if grouping CoMetGene trails by
        reactions)
    """
    rec_limit = 3000  # fails for the default of 1000
    if sys.getrecursionlimit() < rec_limit:
        sys.setrecursionlimit(rec_limit)

    print 'Initialization'

    r_numbers = init_r_numbers(kgml_dir, 1, 4 if by_r else 3)
    reaction_sets = init_reaction_sets(results, 2, 4 if by_r else 3)
    species_dict = init_species_dict(kgml_dir, 3, 4 if by_r else 3)

    if by_r:
        print '\t(4/4) Identifying species reaction sets ...',
        species_reaction_sets = \
            get_species_reaction_sets(reaction_sets, species)
        print 'done'
        return [r_numbers, reaction_sets, species_reaction_sets, species_dict]

    return [r_numbers, reaction_sets, species_dict]


def init_r_numbers(kgml_dir, step, total_steps):
    """Returns a RNumbers object representing associations in KGML files between
    R numbers, species, genes and pathways, for every species in kgml_dir.

    The RNumbers object is retrieved from a pickle file if it exists.

    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :param step: current step during initialization
    :param total_steps: total number of initialization steps
    :return: RNumbers object
    """
    print '\t(%d/%d) Parsing KGML files ...' % (step, total_steps),
    if not os.path.exists(PICKLE_RN):
        r_numbers = get_r_numbers(kgml_dir)
        pickle(PICKLE_RN, r_numbers)
    else:
        r_numbers = unpickle(PICKLE_RN)
        for r_number in r_numbers.r_numbers:
            for organism in r_number.organisms:
                if organism.name not in os.listdir(kgml_dir):
                    raise TrailGroupingError
    print 'done'

    return r_numbers


def get_r_numbers(kgml_dir):
    """Creates and returns a RNumbers object built by parsing associations in
    KGML files between R numbers, species, genes and pathways, for every species
    in kgml_dir.

    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :return: RNumbers object
    """
    r_numbers = RNumbers()

    try:
        for org in os.listdir(kgml_dir):
            org_kgml_dir = os.path.join(kgml_dir, org)
            maps = [f for f in os.listdir(org_kgml_dir) if f.endswith('.kgml')]
            for pathway in maps:
                pw_id = pathway.split('.')[0].split('_')[1].replace(org, '')
                kgml = os.path.join(org_kgml_dir, pathway)
                reactions, compounds = parse_kgml(kgml)
                for r_id in reactions:
                    for rn in reactions[r_id]['reaction']:
                        r_numbers.add(rn.replace('rn:', ''),
                                      org, reactions[r_id]['enzyme'], pw_id)
    except IOError:
        sys.stderr.write('Error reading KGML files.\n')
        exit(1)

    return r_numbers


def init_reaction_sets(results, step, total_steps):
    """Returns a ReactionSets object storing CoMetGeNe trail as CoMetGeNe
    reaction sets.

    The ReactionSets objects is retrieved from a pickle file if it exists.

    :param results: parsed CoMetGeNe results
    :param step: current step during initialization
    :param total_steps: total number of initialization steps
    :return: ReactionSets object
    """
    print '\t(%d/%d) Generating reaction sets ...' % (step, total_steps),
    if not os.path.exists(PICKLE_RS):
        reaction_sets = get_reaction_sets(results)
        pickle(PICKLE_RS, reaction_sets)
    else:
        reaction_sets = unpickle(PICKLE_RS)
    print 'done'

    return reaction_sets


def get_reaction_sets(results):
    """Converts CoMetGeNe trails to CoMetGeNe reaction sets, storing superset
    and subset information.

    :param results: parsed CoMetGeNe results
    :return: ReactionSets object
    """
    reaction_sets = ReactionSets()

    for trail in results.trails:
        reaction_sets.add_set(trail)

    reaction_sets.add_supersets_and_subsets()

    for reaction_set in reaction_sets.sets:
        reaction_set.add_organisms()

    return reaction_sets


def init_species_dict(kgml_dir, step, total_steps):
    """Returns a dictionary storing genomic information for every species in
    kgml_dir.

    The dict is retrieved from a pickle file if it exists.

    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :param step: current step during initialization
    :param total_steps: total number of initialization steps
    :return: dict storing genome information for every species in the data set
        with the exception of the reference species
    """
    print '\t(%d/%d) Building genomes for all species ...' % (
        step, total_steps),
    if not os.path.exists(PICKLE_GEN):
        species_dict = get_species_dict(kgml_dir)
        pickle(PICKLE_GEN, species_dict)
    else:
        species_dict = unpickle(PICKLE_GEN)
        for org in os.listdir(kgml_dir):
            if org not in species_dict:
                raise TrailGroupingError
    print 'done'

    return species_dict


def get_species_dict(kgml_dir):
    """Creates and returns a dictionary storing its genomic information, for
    every species in kgml_dir.

    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :return: dict storing genome information for every species in the data set
        with the exception of the reference species
    """
    species_dict = dict()

    for organism in os.listdir(kgml_dir):
        species_dict[organism] = get_species_object(organism)

    return species_dict


def get_species_object(species, species_genes=None):
    """Creates and returns a Species object with the appropriate Chromosome
    list and gene contents.

    :param species: KEGG organism code
    :param species_genes: genes of 'species' involved in trails common to
        'species' and at least another species; if None, all genes for 'species'
        are used
    :return: a Species object with chromosomes and gene lists for every
        chromosomal strand
    """
    org = Species(species)
    genomes = unpickle(PICKLE_GENOME)
    chromosomes = get_species_chromosomes(genomes, species)

    if species_genes is None and species not in genomes:
        raise TrailGroupingError

    genes = genomes[species] if species_genes is None else species_genes.keys()

    for chromosome in chromosomes:
        current_chr = Chromosome(chromosome)
        org.add_chromosome(current_chr)

        genes_plus = get_genes_on_strand(
            genomes, species, chromosome, genes, True)
        genes_minus = get_genes_on_strand(
            genomes, species, chromosome, genes, False)

        genes_plus = order_genes(genes_plus, species, genomes)
        genes_minus = order_genes(genes_minus, species, genomes)

        for gene in genes_plus:
            current_chr.add_gene_plus(gene)
        for gene in genes_minus:
            current_chr.add_gene_minus(gene)

    return org


def get_species_chromosomes(genomes, species):
    """Returns a list of the species' chromosomes.

    :param genomes: dict of dicts storing gene information for every gene of
        every species in the dict, namely the name of the chromosome on which
        the gene is located, the strand on the chromosome, as well as the
        position of the gene on the chromosome (in nucleotides)
    :param species: KEGG organism code
    :return: list of chromosomes for species
    """
    if species not in genomes:
        raise TrailGroupingError
    chromosomes = set()
    for gene in genomes[species]:
        chromosomes.add(genomes[species][gene]['chr'])
    return list(chromosomes)


def get_genes_on_strand(genomes, species, chromosome, genes, is_plus):
    """Returns the list of genes of the query species on a specific chromosome
    and strand.

    :param genomes: dict of dicts storing gene information for every gene of
        every species in the dict, namely the name of the chromosome on which
        the gene is located, the strand on the chromosome, as well as the
        position of the gene on the chromosome (in nucleotides)
    :param species: KEGG organism code
    :param chromosome: name of chromosome
    :param genes: list of genes of 'species'
    :param is_plus: True for plus strand, False for minus strand
    :return: genes of 'species' among 'genes' that are found on the specified
        strand of the given chromosome
    """
    if species not in genomes:
        raise TrailGroupingError

    genes_chr = [gene for gene in genes
                 if genomes[species][gene]['chr'] == chromosome]
    if is_plus:
        return sorted([gene for gene in genes_chr
                       if genomes[species][gene]['fwd']])
    return sorted([gene for gene in genes_chr
                   if not genomes[species][gene]['fwd']])


def order_genes(genes, species, genomes):
    """Returns a list of ordered genes for the specified species, on a given
    strand of a given chromosome.

    :param genes: list of genes for species (on a given strand of a given
        chromosome)
    :param species: KEGG organism code
    :param genomes: dict of dicts storing gene information for every gene of
        every species in the dict, namely the name of the chromosome on which
        the gene is located, the strand on the chromosome, as well as the
        position of the gene on the chromosome (in nucleotides)
    :return: order of genes on the chromosome strand
    """
    if species not in genomes:
        raise TrailGroupingError

    order = dict()
    for gene in genes:
        order[genomes[species][gene]['pos'][0][0]] = gene

    positions = sorted(order.keys())
    ordered_genes = list()
    for position in positions:
        ordered_genes.append(order[position])

    return ordered_genes


def get_species_reaction_sets(reaction_sets, species):
    """Returns species-specific reaction sets.

    :param reaction_sets: ReactionSets object representing CoMetGeNe reaction
        sets obtained from CoMetGeNe trails
    :param species: query species (KEGG organism code)
    :return: list of species-specific CoMetGeNe reaction sets
    """
    species_reaction_sets = list()

    for reaction_set in reaction_sets.sets:
        if species in reaction_set.organisms:
            maximum_set = True
            for superset in reaction_set.supersets:
                if species in superset.organisms:
                    maximum_set = False
                    break
            if len(reaction_set.supersets) == 0 or maximum_set:
                species_reaction_sets.append(reaction_set)

    return species_reaction_sets
