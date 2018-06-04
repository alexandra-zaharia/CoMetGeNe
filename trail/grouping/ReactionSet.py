"""
These classes are needed for trail grouping, in order to handle CoMetGeNe trails
as reactions sets.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

from trail.utils import split_entities


class ReactionSet(object):
    """Represents a CoMetGeNe trail as a reaction set."""
    def __init__(self, reactions):
        """Creates a ReactionSet object with the specified reactions, and empty
        sets for trails, supersets, subsets, and organisms.

        :param reactions: frozenset of KEGG R numbers
        """
        self.reactions = reactions
        self.trails = set()  # trails having this corresponding ReactionSet
        self.supersets = set()  # include this ReactionSet
        self.subsets = set()  # included in this ReactionSet
        self.organisms = set()  # set of species names having this ReactionSet

    def add_trail(self, trail):
        """Adds the specified trail to the set of trails for this ReactionSet.

        :param trail: Trail object to be added to the set of trails
        """
        self.trails.add(trail)

    def add_superset(self, superset):
        """Adds the specified ReactionSet as superset of this reaction set.

        :param superset: ReactionSet object that includes the reactions of
            this ReactionSet
        """
        self.supersets.add(superset)

    def add_subset(self, subset):
        """Adds the specified ReactionSet as subset of this reaction set.

        :param subset: ReactionSet object whose reactions are included in the
            reactions of this ReactionSet
        """
        self.subsets.add(subset)

    def add_organisms(self):
        """Adds the corresponding organisms to this ReactionSet, i.e. species
        for which the trails corresponding to this ReactionSet occur."""
        for trail in self.trails:
            for organism in trail.organisms:
                self.organisms.add(organism.org)

    def find(self, reactions):
        """Returns this ReactionSet object if it has the specified reactions,
        or None if it is not the case.

        :param reactions: frozenset of KEGG R numbers
        :return: this ReactionSet if its reactions are the ones specified, or
            None otherwise
        """
        return self if self.reactions == reactions else None


class ReactionSets(object):
    """Handles a collection of ReactionSet objects."""
    def __init__(self):
        """Creates ReactionSets object having an empty list of reaction sets."""
        self.sets = list()

    def find(self, reactions):
        """Returns the ReactionSet object having the specified reactions, or
        None if no such reaction set exists in the list of reaction sets.

        :param reactions: frozenset of KEGG R numbers
        :return: the ReactionSet object having the specified reactions, or None
            if no such reaction set is found
        """
        for reaction_set in self.sets:
            rs = reaction_set.find(reactions)
            if rs is not None:
                return rs
        return None

    def add_set(self, trail):
        """Adds a reaction set from the specified Trail object and adds the
        Trail object to the set of trails of the reaction set.

        :param trail: Trail object whose set of reactions becomes a new reaction
            set if no ReactionSet object with the given reactions exists
        """
        reactions = frozenset(trail.reactions)
        rs = self.find(reactions)
        if rs is None:
            rs = ReactionSet(reactions)
            self.sets.append(rs)
        rs.add_trail(trail)

    def add_supersets_and_subsets(self):
        """Determines which reaction sets are supersets for other reaction sets,
        and updates the supersets and subsets for each ReactionSet object
        accordingly.
        """
        for set1 in self.sets:
            for set2 in self.sets:
                if set1 != set2:
                    r1 = set1.reactions
                    r2 = set2.reactions
                    if r1.issubset(r2) and len(r1) < len(r2):
                        set1.add_superset(set2)
                        set2.add_subset(set1)


class SpeciesReactionSets(object):
    """Handles a collection of reaction sets occurring in a given species, as
    well as the genes of the query species involved in the reaction sets."""
    def __init__(self, reaction_sets, species):
        """Creates a SpeciesReactionSets object, with the specified reaction
        sets for a given species.

        A dict (species_genes) is also created. Its keys are genes of the given
        species involved in reaction sets common to 'species' and at least
        another species. The associated values are lists of ReactionSet objects,
        designating reaction sets in which every gene in the dict is involved.

        :param reaction_sets: ReactionSets object representing the collection of
            reaction sets for species
        :param species: KEGG organism code
        """
        self.reaction_sets = reaction_sets
        self.species = species
        self.species_genes = dict()
        self.determine_species_genes()

    def _get_genes(self, trail):
        """Returns a set of genes involved in the given Trail object for the
        query species. If the trail does not occur for the query species, an
        empty set is returned.

        :param trail: Trail object
        :return: set of genes for the query species involved in the given trail,
            or an empty set if the trail does not occur in the query species
        """
        genes_in_trail = set()

        organism = trail.find_by_o_id(self.species)
        if organism is not None:
            for pathway in organism.pathways:
                genes_in_trail = genes_in_trail.union(
                    split_entities(pathway.g_ids))

        return genes_in_trail

    def determine_species_genes(self):
        """Determines the genes in the query species involved in reaction sets
        that are common to species and at least one other organism. The genes
        are keys in the species_genes dict, and the associated values are lists
        of ReactionSet objects representing reaction sets where every gene in
        the dict is involved.

        The 'species_genes' member is modified.
        """
        for reaction_set in self.reaction_sets.sets:
            if self.species in reaction_set.organisms \
                    and len(reaction_set.organisms) > 1:
                for trail in reaction_set.trails:
                    for gene in self._get_genes(trail):
                        if gene not in self.species_genes:
                            self.species_genes[gene] = list()
                        self.species_genes[gene].append(reaction_set)
