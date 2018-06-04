"""
These classes are needed by the CoMetGeNe parser for trail grouping.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""


class ReactionID(object):
    """Handles internal KGML reaction IDs and associates maximum gap parameters.
    """
    def __init__(self, r_id):
        """Creates a ReactionID object with the specified reaction ID and
        maximum gap parameters set to None.

        :param r_id: list of internal KGML reaction IDs
        """
        self.r_id = r_id
        self.max_dG = None
        self.max_dD = None

    def set_max_deltas(self, max_dG, max_dD):
        """Sets maximum gap parameters to the specified values.

        :param max_dG: maximum CoMetGeNe gap parameter for skipping genes
        :param max_dD: maximum CoMetGeNe gap parameter for skipping reactions
        """
        self.max_dG = max_dG
        self.max_dD = max_dD


class Pathway(object):
    """Represents a pathway in which a CoMetGeNe trail occurs by storing the
    KEGG pathway map ID, the list of reaction IDs in the trail, and the list of
    gene names involved in reactions of the trail.
    """
    def __init__(self, p_id, r_ids, g_ids):
        """Creates a Pathway object with the specified ID, reactions, and genes.

        When CoMetGeNe results are parsed, a trail is identified as an ordered
        list of R numbers. For each trail it is known for which species it
        occurs, and in which pathway(s) it occurs. When the internal KGML
        reaction IDs are associated to a given list of R numbers, there can be
        more than one correspondence between R numbers and internal KGML
        reaction IDs. Hence, ir_ids is a list of ReactionID objects instead of a
        single ReactionID object.

        :param p_id: KEGG pathway map ID
        :param r_ids: list of internal KGML reaction IDs associated to a trail
        :param g_ids: list of gene names involved in the reactions in r_ids
        """
        self.p_id = p_id
        self.ir_ids = list()  # list of KGML internal reaction IDs
        self.ir_ids.append(ReactionID(r_ids))
        self.g_ids = g_ids

    def add_r_ids(self, r_ids):
        """If the specified internal KGML reaction IDs are not already found for
        this Pathway object, they are added to the list of ReactionID objects.

        The r_ids member might be modified.

        :param r_ids: list of internal KGML reaction IDs associated to a trail
        """
        if self.find_by_r_id(r_ids) is None:
            self.ir_ids.append(ReactionID(r_ids))

    def find_by_r_id(self, r_id):
        """Returns the ReactionID object corresponding to the given list of
        internal KGML reaction IDs if it exists.

        :param r_id: list of internal KGML reaction IDs associated to a trail
        :return: ReactionID object with the same list of internal KGML reaction
            IDs, or None if no such object exists in this object's ir_ids member
        """
        for ir in self.ir_ids:
            if ir.r_id == r_id:
                return ir
        return None


class Organism(object):
    """Represents a species in which a CoMetGeNe trail occurs."""
    def __init__(self, org):
        """Creates an Organism object with the specified name and with an empty
        list of Pathway objects.

        :param org: name for this Organism
        """
        self.org = org
        self.pathways = list()

    def find_by_p_id(self, p_id):
        """Returns the Pathway object having the specified KEGG pathway map ID
        if it exists in this object's 'pathway' member.

        :param p_id: KEGG pathway map ID
        :return: Pathway object with the same map ID as p_id if it exists in
            this object's 'pathways' member, or None if no such object exists
        """
        for pathway in self.pathways:
            if pathway.p_id == p_id:
                return pathway
        return None

    def add_pathway(self, p_id, r_ids, g_ids):
        """Adds internal KGML reaction IDs and gene names corresponding to these
        reactions to the Pathway object whose KEGG pathway map ID is p_id. If no
        such Pathway object exists in this object's 'pathways' member, creates
        it and adds it to the list.

        The 'pathways' member is modified.

        :param p_id: KEGG pathway map ID
        :param r_ids: list of internal KGML reaction IDs associated to a trail
        :param g_ids: list of gene names involved in the reactions in r_ids
        """
        pathway = self.find_by_p_id(p_id)
        if pathway is None:
            pathway = Pathway(p_id, r_ids, g_ids)
            self.pathways.append(pathway)
        else:
            pathway.add_r_ids(r_ids)


class Trail(object):
    """Represents a CoMetGeNe trail.

    A Trail object is identified by the list of R numbers in the trail. EC
    numbers are additionally stored for the trail, as well as a list of species
    in which the trail occurs. For every species, every pathway in which the
    trail occurs is designated by its KEGG pathway map ID. The internal KGML
    reaction IDs for the given trail of R numbers, species, and pathway are
    equally stored. The case when different internal KGML reaction IDs
    correspond to a given trail of R numbers is handled.
    """
    def __init__(self, reactions, ec_numbers):
        """Creates a Trail object with the given list of R numbers and
        associated EC numbers, and an empty list of Organism objects.

        :param reactions: list of R numbers corresponding to the CoMetGeNe trail
        :param ec_numbers: list of EC numbers corresponding to the CoMetGeNe
            trail
        """
        self.reactions = reactions
        self.ec_numbers = ec_numbers
        self.organisms = list()

    def find_by_o_id(self, org):
        """Returns the Organism object having the specified name if it exists in
        this object's 'organisms' member.

        :param org: species name
        :return: Organism object whose name is designated by org if it exists in
            this object's 'organisms' member, or None if no such object exists
        """
        for organism in self.organisms:
            if organism.org == org:
                return organism
        return None

    def add_organism(self, org):
        """Creates and adds an Organism object whose name is designated by org
        to this object's 'organisms' member, if no such object is present in the
        list.

        The 'organisms' member is modified.

        :param org: species name
        """
        organism = self.find_by_o_id(org)
        if organism is None:
            organism = Organism(org)
            self.organisms.append(organism)

    def add_pathway(self, org, p_id, r_ids, g_ids):
        """Adds a Pathway object to the designated Organism.

        :param org: species name
        :param p_id: KEGG map ID of a pathway in which this Trail occurs for org
        :param r_ids: list of internal KGML reaction IDs associated to this
            Trail object
        :param g_ids: list of gene names involved in the reactions in r_ids for
            this Trail object
        """
        organism = self.find_by_o_id(org)
        assert organism is not None
        organism.add_pathway(p_id, r_ids, g_ids)


class ReactionTrails(object):
    """Represents CoMetGeNe trails (possibly for several species)."""
    def __init__(self):
        """Creates a ReactionTrails object with an empty list of Trail objects.
        """
        self.trails = list()

    def add_trail(self, reactions, ec_numbers):
        """Creates and adds a Trail object with the specified reactions and EC
        numbers to this object's 'trails' member if no such Trail object exists.
        Regardless of whether the Trail object exists initially, it is returned.

        The 'trails' member might be modified.

        :param reactions: list of R numbers corresponding to a CoMetGeNe trail
        :param ec_numbers: list of EC numbers corresponding to a CoMetGeNe trail
        :return: the existing or newly created Trail object
        """
        trail = self.find_by_r_id(reactions)
        if trail is None:
            trail = Trail(reactions, ec_numbers)
            self.trails.append(trail)
        return trail

    def find_by_r_id(self, reactions):
        """Returns the Trail object having the associated specified reactions if
        it exists in this object's 'trails' member.

        :param reactions: list of R numbers corresponding to a CoMetGeNe trail
        :return: Trail object with the R numbers specified by 'reactions' if it
            exists in this object's 'trails' member, or None if no such object
            exists
        """
        for trail in self.trails:
            if trail.reactions == reactions:
                return trail
        return None
