"""
These classes are needed for trail grouping, in order to represent genomic
information for a given species.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import numpy as np


class GeneMatrix(object):
    """Creates 2D boolean arrays."""
    def __init__(self, nb_genes, nb_org):
        """Creates a 2D boolean array of the specified size.

        Used to represent the T_S^g table. For more details, see Trail grouping
        by genes in the methods section of Zaharia et al. (2018).

        :param nb_genes: number of rows for the array
        :param nb_org: number of columns for the array
        """
        self.matrix = np.zeros((nb_genes, nb_org), dtype=bool)


class Chromosome(object):
    """Represents a chromosome of a given species."""
    def __init__(self, name):
        """Creates a new Chromosome object with the specified name and two empty
        gene lists for the two chromosomal strands.

        :param name: name for this chromosome
        """
        self.name = name
        self.plus = list()
        self.minus = list()
        self.matrix_fwd = None  # T_S^g for the + strand of this Chromosome
        self.matrix_rev = None  # T_S^g for the - strand of this Chromosome

    def init_matrix_fwd(self, nb_genes, nb_org):
        """Initializes the matrix for the plus strand of this Chromosome.

        :param nb_genes: number of genes on the plus strand
        :param nb_org: number of species
        """
        self.matrix_fwd = GeneMatrix(nb_genes, nb_org)

    def init_matrix_rev(self, nb_genes, nb_org):
        """Initializes the matrix for the minus strand of this Chromosome.

        :param nb_genes: number of genes on the minus strand
        :param nb_org: number of species
        """
        self.matrix_rev = GeneMatrix(nb_genes, nb_org)

    def get_matrix_fwd(self):
        """Returns the matrix for the plus strand of this chromosome.

        :return: GeneMatrix object for the plus strand
        """
        return None if self.matrix_fwd is None else self.matrix_fwd.matrix

    def get_matrix_rev(self):
        """Returns the matrix for the minus strand of this chromosome.

        :return: GeneMatrix object for the minus strand
        """
        return None if self.matrix_rev is None else self.matrix_rev.matrix

    def add_gene_plus(self, gene):
        """Adds the specified gene to the plus strand of this chromosome.

        :param gene: gene name to add to the plus strand
        """
        self.plus.append(gene)

    def add_gene_minus(self, gene):
        """Adds the specified gene to the minus strand of this chromosome.

        :param gene: gene name to add to the minus strand
        """
        self.minus.append(gene)


class Species(object):
    """Represents genomic information for a given species."""
    def __init__(self, name):
        """Creates a new Species object with the specified name and an empty
        chromosome list.

        :param name: name for this species
        """
        self.name = name
        self.chromosomes = list()

    def add_chromosome(self, chromosome):
        """Adds a Chromosome object to this object's 'chromosomes' member.

        :param chromosome: Chromosome object
        """
        self.chromosomes.append(chromosome)

    def find(self, gene):
        """Returns the chromosome containing the specified gene (if it exists).

        :param gene: query gene
        :return: Chromosome object containing the query gene or None if no such
            chromosome exists
        """
        for chromosome in self.chromosomes:
            if gene in chromosome.plus or gene in chromosome.minus:
                return chromosome
        return None

    def find_chromosome(self, chr_name):
        """Returns the chromosome of this species with the specified name (if it
        exists).

        :param chr_name: name of chromosome to search for this species
        :return: Chromosome object with the same name as chr_name or None if no
            such chromosome exits
        """
        for chromosome in self.chromosomes:
            if chromosome.name == chr_name:
                return chromosome
        return None

    def is_plus(self, gene, chromosome):
        """Checks whether the query gene is present on the plus or minus strand
        on the specified chromosome. Chromosome must exist for this species.

        :param gene: query gene
        :param chromosome: Chromosome on which the gene must be found
        :return: True if gene is on the plus strand, False otherwise
        """
        assert chromosome in self.chromosomes
        return gene in chromosome.plus
