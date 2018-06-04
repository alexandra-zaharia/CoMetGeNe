from trail.utils import split_entities


class NeighborhoodResolver(object):
    """Resolves neighboring reactions (used in trail grouping by reactions)."""
    DELTA = 4  # maximum number of skipped genes + 1

    ABSENT = 'o'  # symbol for absent reaction
    NON_CONTIG = '.'  # symbol for reaction with non contiguous gene(s)
    CONTIG = 'x'  # symbol for reaction with contiguous gene(s)

    def __init__(self, species, reaction_set, r_numbers):
        """Initializes this NeighborhoodResolver instance with the necessary
        information in order to determine a maximal reaction neighborhood for
        the target species.

        :param species: Species object representing a target species other than
            the reference species
        :param reaction_set: ReactionSet object representing a reaction set
            for the target species for which a maximal reaction neighborhood
            will be determined
        :param r_numbers: RNumbers object associating R numbers, species, genes,
            and pathways
        """
        self.species = species
        self.reaction_set = reaction_set
        self.r_numbers = r_numbers

        self.genes_to_reactions = dict()
        self.gene_list = self._get_gene_list()

        self.absent_reactions = set()
        self._get_absent_reactions()

        self.indexes = self._get_indexes()

        self.neighborhood_candidates = dict()
        self._get_neighborhood_candidates()

        self.best = self.best_candidate()
        if self.best is not None:
            self.bc_reactions = self._get_candidate_reactions(self.best)

    def _get_gene_list(self):
        """Returns the list of genes for the target species involved in
        reactions in this object's reaction_set.

        :return: list of genes of the target species involved in this object's
            reaction set
        """
        gene_list = list()

        for reaction in split_entities(self.reaction_set.reactions):
            r_number = self.r_numbers.find(reaction)
            assert r_number is not None
            organism = r_number.find(self.species.name)
            if organism is not None:  # else reaction absent in species
                gene_list.append(organism.genes)
                for gene in organism.genes:
                    if gene not in self.genes_to_reactions:
                        self.genes_to_reactions[gene] = list()
                    self.genes_to_reactions[gene].append(reaction)

        return gene_list

    def _get_absent_reactions(self):
        """Adds the reactions which are absent for the target species to this
        object's absent_reactions member."""
        for reaction in self.reaction_set.reactions:
            for r_number in split_entities([reaction]):
                organism = self.r_numbers.find(r_number).find(self.species.name)
                if organism is None:
                    self.absent_reactions.add(reaction)

    def _get_indexes(self):
        """Determines the indexes of genes involved in the reaction set, for
        every chromosome and every strand of the target species.

        :return: dict of dicts with chromosome names as keys, strands as
            sub-keys, gene indexes on the given chromosome and strand as
            sub-sub-keys, and finally gene names as values
        """
        indexes = dict()

        for gene_set in self.gene_list:
            for gene in gene_set:
                chromosome = self.species.find(gene)
                assert chromosome is not None
                is_plus = self.species.is_plus(gene, chromosome)
                index = chromosome.plus.index(gene) \
                    if is_plus else chromosome.minus.index(gene)
                if chromosome.name not in indexes:
                    indexes[chromosome.name] = dict()
                strand = '+' if is_plus else '-'
                if strand not in indexes[chromosome.name]:
                    indexes[chromosome.name][strand] = dict()
                indexes[chromosome.name][strand][index] = gene

        return indexes

    def _empty_valid_helper(self, indexes, i, valid):
        """Adds indexes i and i+1 to the list of valid indexes if the difference
        between the two indexes is at most DELTA.

        The 'valid' list is modified.

        :param indexes: list of sorted indexes
        :param i: index in the list of sorted indexes
        :param valid: list of valid indexes in terms of neighborhood capability
        """
        if i < len(indexes) - 1 and indexes[i + 1] - indexes[i] <= self.DELTA:
            valid.append(indexes[i])
            valid.append(indexes[i + 1])

    def _full_valid_helper(self, chromosome, strand, valid):
        """Adds the current valid list of indexes (in terms of neighborhood
        capability) to the neighborhood_candidates member of this object, for
        the given chromosome and strand.

        :param chromosome: chromosome name for the possible neighborhood
        :param strand: strand symbol ('+' or '-') for the possible neighborhood
        :param valid: list of valid indexes in terms of neighborhood capability
        """
        self.neighborhood_candidates[frozenset(valid)] = dict()
        self.neighborhood_candidates[frozenset(valid)]['chr'] = chromosome
        self.neighborhood_candidates[frozenset(valid)]['strand'] = strand

    def _get_neighborhood_candidates(self):
        """Identifies neighborhood candidates, i.e. groups of gene indexes for
        the same chromosome and the same strand such that they differ by at most
        DELTA indexes.

        This object's neighborhood_candidates member is updated.
        """
        for chromosome in self.indexes:
            for strand in self.indexes[chromosome]:
                sorted_idx = sorted(
                    [index for index in self.indexes[chromosome][strand]])

                valid = list()

                for i in range(len(sorted_idx)):
                    if len(valid) == 0:
                        self._empty_valid_helper(sorted_idx, i, valid)
                    elif sorted_idx[i] != valid[-1]:
                        if sorted_idx[i] - valid[-1] <= self.DELTA:
                            valid.append(sorted_idx[i])
                        else:
                            self._full_valid_helper(chromosome, strand, valid)
                            del valid[:]
                            self._empty_valid_helper(sorted_idx, i, valid)

                if len(valid) > 0:
                    self._full_valid_helper(chromosome, strand, valid)

    def _get_candidate_reactions(self, candidate):
        """Returns the set of reactions corresponding to the indexes in
        'candidate' from this object's 'indexes' member.

        :param candidate: frozenset of gene indexes
        :return: set of R numbers of corresponding reactions for the 'candidate'
            genes designated by this object's 'indexes' member
        """
        reactions = set()

        chromosome = self.neighborhood_candidates[candidate]['chr']
        strand = self.neighborhood_candidates[candidate]['strand']

        for index in sorted(candidate):
            gene = self.indexes[chromosome][strand][index]
            for reaction in self.genes_to_reactions[gene]:
                reactions.add(reaction)

        return reactions

    def best_candidate(self):
        """Chooses the best neighborhood candidate, i.e. the one that maximizes
        the number of reactions from the original reaction set in which its
        genes are involved.

        :return: frozenset of indexes designating the best neighborhood
            candidate among candidates in this object's neighborhood_candidates
            member
        """
        best_candidate = None
        number_of_reactions = 0

        for candidate in self.neighborhood_candidates:
            reactions = self._get_candidate_reactions(candidate)
            if number_of_reactions < len(reactions):
                number_of_reactions = len(reactions)
                best_candidate = candidate

        return best_candidate

    def state(self, reaction):
        """Determines whether a given reaction is present, in which case the
        distinction is made between reactions catalyzed by products of a gene
        neighboring another gene involved in other reactions of the current
        reaction trail on the one hand, and reactions catalyzed by products of
        genes which are not neighbors of other genes involved in the same
        reaction trail, on the other hand.

        :param reaction: KEGG R number
        :return: symbol describing whether the reaction is present or absent
        """
        if reaction in self.absent_reactions:
            return self.ABSENT

        if self.best is None:
            return self.NON_CONTIG

        for r_number in split_entities([reaction]):
            if r_number not in self.bc_reactions:
                return self.NON_CONTIG

        return self.CONTIG


class NeighborhoodResolverGenes(NeighborhoodResolver):
    """Resolves neighboring genes (used in trail grouping by genes)."""
    def __init__(self, reference, species, reaction_set, r_numbers):
        """Initializes this NeighborhoodResolver instance with the necessary
        information in order to determine a maximal gene neighborhood for the
        target species.

        :param reference: Species object representing the reference species
        :param species: Species object representing a target species other than
            the reference species
        :param reaction_set: ReactionSet object representing a reaction set
            for the target species for which a maximal reaction neighborhood
            will be determined
        :param r_numbers: RNumbers object associating R numbers, species, genes,
            and pathways
        """
        self.reference = reference
        super(NeighborhoodResolverGenes, self)\
            .__init__(species, reaction_set, r_numbers)

    def best_candidate(self):
        """Chooses the best neighborhood candidate, i.e. the one that maximizes
        the number of genes of the reference species having neighboring
        functionally similar genes in the target species being involved in
        reactions from the original reaction set.

        :return: frozenset of indexes designating the best neighborhood
            candidate among candidates in this object's neighborhood_candidates
            member
        """
        best_candidate = None
        number_of_genes = 0

        for candidate in self.neighborhood_candidates:
            reactions = self._get_candidate_reactions(candidate)
            gene_set = set()
            for reaction in reactions:
                r_number = self.r_numbers.find(reaction)
                assert r_number is not None
                ref_org = r_number.find(self.reference.name)
                assert ref_org is not None
                for gene in ref_org.genes:
                    gene_set.add(gene)

            if number_of_genes < len(gene_set):
                number_of_genes = len(gene_set)
                best_candidate = candidate

        return best_candidate
