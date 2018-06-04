"""
These methods are used for grouping CoMetGeNe trails by reactions.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import re

import sys

from NeighborhoodResolver import NeighborhoodResolverGenes
from ReactionSet import SpeciesReactionSets, ReactionSet
from ..definitions import MAX_DELTA, PICKLE_GENOME
from trail.utils import unpickle
from trail_grouping_utils import get_species_object, order_genes, \
    get_species_chromosomes, get_genes_on_strand, initialization
from phylogeny import get_ordered_organisms


def fill_gene_matrix(results, species, kgml_dir):
    """Returns a Species object for the reference species with required genomic
    information, as well as T_S^g tables for every chromosome.

    :param results: parsed CoMetGeNe results
    :param species: reference species (KEGG organism code)
    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :return: a Species object with filled chromosomes and gene lists for every
        chromosomal strand, as well as T_S^g tables for every chromosome (for
        more details, see Trail grouping by genes in the methods section of
        Zaharia et al., 2018)
    """
    r_numbers, reaction_sets, species_dict = initialization(
        results, species, kgml_dir, False)

    species_sets = SpeciesReactionSets(reaction_sets, species)
    species_genes = species_sets.species_genes

    organisms = get_ordered_organisms(species, kgml_dir)
    organism = get_species_object(species, species_genes=species_genes)
    remove_non_neighboring_genes(organism)

    for chromosome in organism.chromosomes:
        chromosome.init_matrix_fwd(len(chromosome.plus), len(organisms))
        chromosome.init_matrix_rev(len(chromosome.minus), len(organisms))
        for j in range(len(organisms)):
            resolve_neighboring(organism, chromosome, True,
                                species_sets, r_numbers,
                                j,
                                species_dict[organisms[j]])
            resolve_neighboring(organism, chromosome, False,
                                species_sets, r_numbers,
                                j,
                                species_dict[organisms[j]])

    return organism


def remove_non_neighboring_genes(organism):
    """Removes genes of the given Species object 'organism', from every
    chromosome and strand, if these genes are 'singletons', i.e. if they are not
    neighbors of other genes present in the Species object.

    Note that singleton genes are almost always present in Chromosome objects of
    Species objects. Recall that a Species object stores genomic information for
    genes involved in CoMetGeNe trails common to the reference species (for
    which a Species object is created) and at least one other species in the
    data set.

    Chromosome objects of the Species object 'organism' are modified.

    :param organism: Species object with chromosomes and gene lists for every
        chromosomal strand
    """
    genomes = unpickle(PICKLE_GENOME)
    genes = genomes[organism.name]
    chromosomes = get_species_chromosomes(genomes, organism.name)

    for chromosome in chromosomes:
        genes_plus = get_genes_on_strand(
            genomes, organism.name, chromosome, genes, True)
        genes_minus = get_genes_on_strand(
            genomes, organism.name, chromosome, genes, False)
        ordered_plus = order_genes(genes_plus, organism.name, genomes)
        ordered_minus = order_genes(genes_minus, organism.name, genomes)

        org_chr = organism.find_chromosome(chromosome)
        assert org_chr is not None

        remove_plus = determine_genes_to_remove(org_chr.plus, ordered_plus)
        remove_minus = determine_genes_to_remove(org_chr.minus, ordered_minus)

        for gene in remove_plus:
            org_chr.plus.remove(gene)
        for gene in remove_minus:
            org_chr.minus.remove(gene)

        org_chr.groups_plus = delimit_gene_groups(org_chr.plus, ordered_plus)
        org_chr.groups_minus = delimit_gene_groups(org_chr.minus, ordered_minus)


def determine_genes_to_remove(org_strand, genome_strand):
    """Determines singleton genes on a chromosomal strand, i.e. genes that are
    not neighbors of other genes involved in CoMetGeNe trails.

    :param org_strand: list of gene identifiers on a given strand of a given
        chromosome, designating genes involved in CoMetGeNe trails
    :param genome_strand: list of genes identifiers on a given strand of a given
        chromosome in genomic order, representing all the genes present on the
        strand (irrespective of their involvement in CoMetGeNe trails)
    :return: list of singleton genes on the strand 'org_strand', i.e. genes
        separated by more than MAX_DELTA genes from the next gene in org_strand
    """
    genes_to_remove = list()

    for i in range(len(org_strand)):
        remove = True

        assert org_strand[i] in genome_strand
        idx_curr = genome_strand.index(org_strand[i])

        if i == 0:
            idx_prev = idx_curr
            if len(org_strand) > 1:
                idx_next = genome_strand.index(org_strand[1])
            else:
                idx_next = idx_curr
        elif i == len(org_strand) - 1:
            idx_next = idx_curr
            if len(org_strand) > 1:
                idx_prev = genome_strand.index(org_strand[i-1])
            else:
                idx_prev = idx_curr
        else:
            idx_prev = genome_strand.index(org_strand[i-1])
            idx_next = genome_strand.index(org_strand[i+1])

        if idx_next != idx_curr and idx_next - idx_curr <= MAX_DELTA + 1:
            remove = False
        elif idx_curr != idx_prev and idx_curr - idx_prev <= MAX_DELTA + 1:
            remove = False

        if remove:
            genes_to_remove.append(org_strand[i])

    return genes_to_remove


def delimit_gene_groups(org_strand, genome_strand):
    """Returns lists of indexes of neighboring genes in org_strand, where the
    indexes designate the positions of neighboring genes in genome_strand.

    Genes on the same strand of a given chromosome are considered neighbors if
    they are separated by at most MAX_DELTA other genes. For example, for the
    following input and MAX_DELTA = 3:

        org_strand        B     D  E              J  K  L     N
        genome_strand  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P
        (index         0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15)

    the method returns

        [[1, 3, 4], [9, 10, 11, 13]]

    because the genes B, D, E and J, K L, N, respectively, are neighbors, since
    they are separated by less than MAX_DELTA other genes. On the contrary,
    genes E and J are not neighbors, because they are separated by MAX_DELTA + 1
    other genes.

    :param org_strand: list of gene identifiers on a given strand of a given
        chromosome, designating genes involved in CoMetGeNe trails
    :param genome_strand: list of genes identifiers on a given strand of a given
        chromosome in genomic order, representing all the genes present on the
        strand (irrespective of their involvement in CoMetGeNe trails)
    :return: list of lists, where each list represents indexes of neighboring
        genes in org_strand on the genome_strand
    """
    gene_groups = list()
    current_group = list()

    for i in range(len(org_strand)):
        idx_curr = genome_strand.index(org_strand[i])
        current_group.append(i)

        idx_next = idx_curr
        if i < len(org_strand) - 1:
            idx_next = genome_strand.index(org_strand[i+1])

        if idx_curr == idx_next or idx_next - idx_curr > MAX_DELTA + 1:
            gene_groups.append(current_group[:])
            del current_group[:]

    return gene_groups


def resolve_neighboring(
        reference,
        chromosome, is_plus,
        species_sets, r_numbers, j, other_species):
    """Resolves neighboring for the j-th species in the data set for genes
    functionally similar to genes of the reference species, by filling column j
    in T_S^r. For more details, see Trail grouping by genes in the methods
    section of Zaharia et al. (2018).

    The GeneMatrix member of the Chromosome object 'chromosome' of the reference
    species is modified for the given strand (matrix_fwd for plus strand,
    matrix_rev for minus strand).

    :param reference: Species object representing the reference species
    :param chromosome: Chromosome object of the reference species
    :param is_plus: True for plus strand of chromosome, False for minus strand
    :param species_sets: SpeciesReactionSets object representing CoMetGeNe
    reaction sets for the reference species, obtained from CoMetGeNe trails
    :param r_numbers: RNumbers object associating R numbers, species, genes and
        pathways, for every species in the data set, as determined from KGML
        files
    :param j: index of column in T_S^g to fill, representing the j-th species
        in the data set sorted by phylogenetic order (or alphabetical order if
        the data set contains species other than those in Table 2 of Zaharia et
        al., 2018 and the user did not yet define a phylogenetic order for the
        new data set)
    :param other_species: Species object representing a target species other
        than the reference species
    """
    groups = chromosome.groups_plus if is_plus else chromosome.groups_minus
    genes_chr = chromosome.plus if is_plus else chromosome.minus
    matrix = chromosome.get_matrix_fwd() if is_plus \
        else chromosome.get_matrix_rev()

    for group in groups:
        genes = [genes_chr[index] for index in group]
        g2r, r2g = associate_genes_and_r_numbers(
            genes, species_sets.species, r_numbers)
        reactions = get_unique_r_numbers(g2r)

        n_resolver = NeighborhoodResolverGenes(
            reference, other_species, ReactionSet(reactions), r_numbers)
        for reaction in reactions:
            if n_resolver.state(reaction) == n_resolver.CONTIG:
                for gene in r2g[reaction]:
                    if gene in genes:
                        matrix[genes_chr.index(gene), j] = True


def associate_genes_and_r_numbers(genes, species, r_numbers):
    """Associates genes and R numbers for the reference species.

    For a given R number, the associated genes are genes of the reference
    species encoding enzymes that catalyze a reaction with the given R number.

    For a given gene, the associated R numbers designate reactions which are
    catalyzed by the product of the given gene.

    :param genes: list of gene identifiers of the reference species
    :param species: KEGG organism code for the reference species
    :param r_numbers: RNumbers object associating R numbers, species, genes and
        pathways, for every species in the data set, as determined from KGML
        files
    :return: two dicts, the first representing R numbers associated to a given
        gene name, the second representing gene names associated to a given R
        number
    """
    g2r = dict()
    r2g = dict()

    for gene in genes:
        g2r[gene] = set()
        for r_number in r_numbers.r_numbers:
            org_obj = r_number.find(species)
            if org_obj is not None and gene in org_obj.genes:
                    reaction = r_number.r_number
                    g2r[gene].add(reaction)
                    if reaction not in r2g:
                        r2g[reaction] = set()
                    r2g[reaction].add(gene)

    return g2r, r2g


def get_unique_r_numbers(g2r):
    """Returns the set of R numbers associated to at least one gene in the g2r
    dict.

    :param g2r: dict with gene names for keys and the associated R numbers for
        values
    :return: frozenset of all R number values in g2r
    """
    reactions = set()

    for gene in g2r:
        for reaction in g2r[gene]:
            reactions.add(reaction)

    return frozenset(reactions)


def print_header(species, organisms, dev_out):
    """Prints table header for the reference species to the specified output
    device.

    :param species: KEGG organism code for the reference species
    :param organisms: list of other species in the data set
    :param dev_out: output device (file handle or stdout)
    """
    sep = ';' if dev_out != sys.stdout else '  '
    header = species + '_gene' + sep + 'chr' + sep + 'str' + sep
    for org in organisms:
        header += org
        if organisms.index(org) < len(organisms) - 1:
            header += ';' if dev_out != sys.stdout is not None else ' '
    dev_out.write(header + '\n')


def end_of_group_reached(i, groups):
    """Determines whether the end of a group of neighboring genes is reached
    when displaying the table associated to trail grouping by genes.

    :param i: index of a gene on a chromosomal strand
    :param groups: list of lists of indexes of neighboring genes on a
        chromosomal strand
    :return: True if i is equal to the last value (index) of a list in 'groups',
        False if i is less than the first value (index) of a list in 'groups'
    """
    for group in groups:
        if i == group[-1]:
            return True
        if i < group[0]:
            return False


def print_species_table(organism, kgml_dir, dev_out):
    """Prints the table corresponding to trail grouping by reactions for the
    reference species to the output device dev_out.

    :param organism: Species object with filled chromosomes and gene lists for
        every chromosomal strand, as well as T_S^g tables for every chromosome
    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :param dev_out: output device (file handle or stdout)
    """
    species = organism.name
    organisms = get_ordered_organisms(species, kgml_dir)
    print_header(species, organisms, dev_out)

    for chromosome in organism.chromosomes:
        plus = chromosome.plus
        minus = chromosome.minus
        for gene_set in plus, minus:
            groups = chromosome.groups_plus if gene_set == plus \
                else chromosome.groups_minus
            for i in range(len(gene_set)):
                gene_line = plus[i] if gene_set == plus else minus[i]
                gene_line = re.sub(species + ':', '', gene_line)
                sep = ';' if dev_out != sys.stdout is not None else '    '
                strand = '+' if gene_set == plus else '-'
                chr_name = re.sub('chromosome', 'chr', chromosome.name)
                gene_line += sep + chr_name + sep + strand
                for j in range(len(organisms)):
                    cell = chromosome.get_matrix_fwd()[i, j] \
                        if gene_set == plus \
                        else chromosome.get_matrix_rev()[i, j]
                    if dev_out != sys.stdout:
                        gene_line += ';x' if cell else ';.'
                    else:
                        gene_line += ' x ' if cell else '   '
                dev_out.write(gene_line + '\n')
                if end_of_group_reached(i, groups):
                    dev_out.write('***\n')

    if dev_out != sys.stdout:
        dev_out.close()
