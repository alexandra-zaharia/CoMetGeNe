"""
These classes are needed for trail grouping, in order to handle associations in
KGML files between R Numbers, species, genes, and pathways.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""


class Organism(object):
    """Stores species-related information for a given R number.

    The information being stored is the species name, the set of genes for this
    species that are involved in the given R number, and the set of pathways
    where this R number is present for the species.
    """
    def __init__(self, name):
        """Creates a new Organism object with the given species name, and empty
        sets for genes and pathways.

        :param name: name for this Organism
        """
        self.name = name
        self.genes = set()
        self.pathways = set()

    def add_gene(self, gene):
        """Adds the given gene to the set of genes of this Organism.

        :param gene: gene to add to the genes set
        """
        self.genes.add(gene)

    def add_pathway(self, pathway):
        """Adds the given pathway to the set of pathways of this Organism.

        :param pathway: pathway to add to the pathways set
        """
        self.pathways.add(pathway)


class RNumber(object):
    """Describes a KEGG R number.

    Each R number has an ID and a set of Organism objects. The Organism objects
    describe which genes in the given species are involved in a particular R
    number and in which pathways the R number is found for the given species.
    """
    def __init__(self, r_number):
        """Creates a new RNumber object with the specified ID and an empty set
        of organisms.

        :param r_number: KEGG R number for this object
        """
        self.r_number = r_number
        self.organisms = set()

    def add_organism(self, org):
        """Adds a new species to this object's set of Organism objects.

        :param org: string for Organism name
        :return: the newly created Organism object
        """
        species = Organism(org)
        self.organisms.add(species)
        return species

    def find(self, org):
        """Returns the Organism object having the specified name, or None if no
        such object exists in the set of organisms for this RNumber.

        :param org: organism name
        :return: the Organism having the specified name, or None if no such
            organism exists in the organisms set
        """
        for organism in self.organisms:
            if organism.name == org:
                return organism
        return None


class RNumbers(object):
    """Stores all associations in KGML files between R Numbers, species, genes,
    and pathways."""
    def __init__(self):
        """Creates a new RNumbers object an empty set of RNumber objects."""
        self.r_numbers = set()

    def find(self, rn_id):
        """Returns the RNumber object having the specified KEGG R number, or
        None if no such object exists in the set of R numbers.

        :param rn_id: KEGG R number
        :return: the RNumber having the specified ID, or None if no such RNumber
            exists in the set
        """
        for r_number in self.r_numbers:
            if r_number.r_number == rn_id:
                return r_number
        return None

    def add(self, rn_id, org, genes, pathway):
        """Creates or updates an existing RNumber with the species where the
        KEGG R number occurs, along with the list of genes involved in the R
        number as well as the pathway in which the R number occurs.

        :param rn_id: KEGG R number
        :param org: species name (KEGG organism code) where the R number appears
        :param genes: list of genes of org involved in the R number
        :param pathway: KEGG pathway map ID where genes of org are involved in
            the R number
        :return: RNumber object (created or updated)
        """
        rn_obj = self.find(rn_id)
        if rn_obj is None:
            rn_obj = RNumber(rn_id)
            self.r_numbers.add(rn_obj)
        org_obj = rn_obj.find(org)
        if org_obj is None:
            org_obj = rn_obj.add_organism(org)
        for gene in genes:
            org_obj.add_gene(gene)
        org_obj.add_pathway(pathway)
        return rn_obj
