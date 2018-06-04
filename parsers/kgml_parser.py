"""
These utilities handle KGML files.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

from lxml import etree


def get_pathway_id(kgml):
    """Returns the KEGG 5-digit ID of a given metabolic pathway in KGML format.

    :param kgml: filename for input metabolic pathway
    :return: string representing the KEGG 5-digit ID of kgml
    """
    tree = etree.parse(kgml)
    root = tree.getroot()
    assert 'number' in root.attrib
    return root.get('number')


def get_pathway_title(kgml):
    """Returns the title of a given metabolic pathway in KGML format.

    :param kgml: filename for input metabolic pathway
    :return: string representing the title of kgml
    """
    tree = etree.parse(kgml)
    root = tree.getroot()
    assert 'title' in root.attrib
    return root.get('title')


def parse_kgml(kgml):
    """Parses the KGML file returns two dictionaries, respectively storing
    information on reactions and compounds in kgml.

    The first dictionary, 'reactions', is actually a dictionary of dictionaries,
    storing reaction identifiers (as they appear in the kgml file) as keys that
    have another dictionary as associated value; this second dictionary (the
    associated value) contains information on enzyme names, as well as R
    numbers, substrate and product identifiers, along with reaction
    reversibility.

    The second dictionary, 'compounds', associates compound names (values) to
    compound identifiers as they appear in kgml.

    :param kgml: filename for input metabolic pathway
    :return: a dict of dicts storing reaction information and a dict storing
        compound information
    """
    tree = etree.parse(kgml)
    root = tree.getroot()

    # 'reactions' is a dictionary of dictionaries: to each key 'r_id' in the
    # dictionary is assigned a new dictionary, with the following entries/keys:
    #     * enzyme     - name(s) of enzyme(s) catalyzing the reaction
    #     * reaction   - KEGG identifier(s) of the reaction designated by r_id
    #     * substrate  - list of identifiers for compounds serving as reaction
    #                    substrates
    #     * product    - list of identifiers for compounds serving as reaction
    #                    products;
    #     * reversible - True if the reaction is reversible, False otherwise
    reactions = dict()

    # 'compounds' is a simple dictionary, associating compound names (values) to
    # their identifiers (keys) as given in the KGML file.
    compounds = dict()

    for elem in root:
        kgml_id = elem.get('id')

        if elem.tag == 'entry':
            if elem.get('type') == 'gene':
                reactions[kgml_id] = dict()

                reactions[kgml_id]['enzyme'] = list()
                names = elem.get('name')
                for name in names.split():
                    reactions[kgml_id]['enzyme'].append(name)

                reactions[kgml_id]['reaction'] = list()
                if 'reaction' in elem.attrib:
                    reaction_ids = elem.get('reaction')
                    for reaction in reaction_ids.split():
                        reactions[kgml_id]['reaction'].append(reaction)

                reactions[kgml_id]['substrate'] = list()
                reactions[kgml_id]['product'] = list()
            elif elem.get('type') == 'compound':
                compounds[kgml_id] = elem.get('name')

        elif elem.tag == 'reaction' and elem.get('id') in reactions:
            reactions[kgml_id]['reversible'] = True \
                if elem.get('type') == 'reversible' else False
            for cpd in elem:
                reactions[kgml_id][cpd.tag].append(cpd.get('id'))

    return reactions, compounds


def add_ec_numbers(reactions, ec_numbers):
    """Adds EC number associations to R numbers in 'reactions'.

    The 'reactions' dict is modified.

    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :param ec_numbers: dict associating R numbers (keys) to EC numbers (values)
    """
    for r_id in reactions:
        reactions[r_id]['ec'] = list()
        for rn in reactions[r_id]['reaction']:
            if rn in ec_numbers:
                for ec in ec_numbers[rn]:
                    reactions[r_id]['ec'].append(ec)
