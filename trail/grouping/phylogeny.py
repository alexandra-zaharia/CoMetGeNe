"""
Utilities for ordering analyzed species in phylogenetic order for trail grouping
or, failing that, alphabetical order.

If trail grouping is executed on other species than those in Table 2 of Zaharia
et al. (2018), the user is warned that the chosen species are displayed in
alphabetical order and is kindly invited to modify methods phylogeny_user and
get_ordered_organisms_user in this file to reflect phylogenetic relations
between species.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import os
import sys


def warning_phylogeny(kgml_dir):
    """Outputs a warning to stderr if kgml_dir contains species not found in
    Table 2 of Zaharia et al. (2018).

    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    """
    if not is_paper_phylogeny_used(kgml_dir):
        sys.stderr.write(
            "\n#############\n"
            "#  WARNING\n"
            "#############\n#\n"
            "# " + kgml_dir + " contains species not found in Table 2 of\n"
            "# Zaharia et al. (2018). In order for trail grouping to contain\n"
            "# species in your data set in phylogenetic order, please modify\n"
            "# methods phylogeny_user and get_ordered_organisms_user in\n"
            "# grouping/phylogeny.py.\n\n"
        )


def phylogeny_user(kgml_dir):
    """Unless modified, this method returns the species under kgml_dir in
    alphabetical order.

    Every subdirectory under kgml_dir is assumed to be the three- or four-letter
    KEGG code for a species. CoMetGeNe users should modify this method and
    get_ordered_organisms_user to return the list of species under study in
    phylogenetic order.

    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :return: list of species under study in alphabetical order
    """
    return sorted(os.listdir(kgml_dir))


def phylogeny_paper():
    """Returns the 50 bacterial species in Table 2 (Zaharia et al., 2018) in
    phylogenetic order.

    :return: list of species under study in phylogenetic order
    """
    proteobacteria = ['eco', 'ype', 'vco', 'spc', 'pae', 'xfa', 'rso', 'nme',
                      'afi', 'ara', 'rrj', 'gsu']
    inter_proteo_terra = ['nde', 'aca', 'din', 'fnu', 'dap', 'tid', 'aae']
    terrabacteria = ['bsu', 'lmo', 'sau', 'lac', 'snd', 'cpe', 'mpn',  # Firm.
                     'syn', 'pma',  # Cyanobacteria
                     'cau',  # Chloroflexi
                     'bbv', 'cgl', 'mtv', 'sco',  # Actinobacteria
                     'dra', 'tth',  # Thermi
                     'fgi']  # Armatimonadetes
    inter_terra_fcb = ['amo', 'tmm', 'cex', 'dth']
    fcb_bacteria = ['fsu', 'gau', 'cph', 'bfr']
    pvc_bacteria = ['rba', 'cpn', 'ote']
    post_pvc = ['bbn', 'emi', 'heo']

    phylogeny = list()
    for bacteria in proteobacteria, inter_proteo_terra, terrabacteria, \
            inter_terra_fcb, fcb_bacteria, pvc_bacteria, post_pvc:
        phylogeny.extend(bacteria)

    return phylogeny


def is_paper_phylogeny_used(kgml_dir):
    """Determines whether the 50 bacterial species in Table 2 (Zaharia et al.,
    2018) should be returned according to the contents of the kgml_dir
    directory.

    If studying other species than those in Table 2, CoMetGeNe users should
    modify phylogeny_user and get_ordered_organisms_user to return the list of
    species under study in phylogenetic order.

    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :return: True if all species from Table 2 (Zaharia et al., 2018) and none
        other are used, False otherwise
    """
    return sorted(phylogeny_paper()) == sorted(os.listdir(kgml_dir))


def get_ordered_organisms(species, kgml_dir):
    """Returns a list of organisms ordered according to their phylogenetic
    distance to the query species.

    :param species: KEGG organism code for the reference species
    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :return: list of organisms, the first being the closest to species from a
        phylogenetic point of view
    """
    if is_paper_phylogeny_used(kgml_dir):
        return get_ordered_organisms_paper(species)
    return get_ordered_organisms_user(species, kgml_dir)


def get_ordered_organisms_user(species, kgml_dir):
    """Unless modified, this method returns the list of organisms under kgml_dir
    in alphabetical order.

    If studying other species than those in Table 2 (Zaharia et al., 2018),
    CoMetGeNe users should modify phylogeny_user and this method to return the
    list of species under study in phylogenetic order.

    :param species: KEGG organism code for the reference species
    :param kgml_dir: directory with a subdirectory for every species that
        CoMetGeNe was ran on, where each subdirectory contains metabolic
        pathways for the species in question in KGML format
    :return: list of organisms, the first being the closest to species from a
        phylogenetic point of view
    """
    phylogeny = phylogeny_user(kgml_dir)

    if species not in phylogeny:
        sys.stderr.write('Undefined tree for species \'' + species + '\'.\n')
        sys.exit(1)

    phylogeny.remove(species)
    return phylogeny


def get_ordered_organisms_paper(species):
    """Returns the list of organisms corresponding to species in Table 2 of
    Zaharia et al. (2018), ordered according to their phylogenetic distance to
    the query species.

    :param species: KEGG organism code for the reference species
    :return: list of organisms, the first being the closest to species from a
        phylogenetic point of view
    """
    phylogeny = phylogeny_paper()

    if species not in phylogeny:
        sys.stderr.write('Undefined tree for species \'' + species + '\'.\n')
        sys.exit(1)

    new_phylogeny = list()
    idx = phylogeny.index(species)

    if species in ['eco', 'ype', 'vco', 'spc', 'pae', 'xfa']:
        new_phylogeny = phylogeny
        new_phylogeny.remove(species)
    if species in ['rso', 'ara']:
        new_phylogeny.extend([phylogeny[idx + 1]])
        new_phylogeny.extend(phylogeny[:idx])
        new_phylogeny.extend(phylogeny[idx + 2:])
    elif species in ['nme', 'rrj']:
        new_phylogeny.extend([phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[:idx - 1])
        new_phylogeny.extend(phylogeny[idx + 1:])
    elif species in ['afi', 'gsu']:
        new_phylogeny.extend(phylogeny[:idx])
        new_phylogeny.extend(phylogeny[idx + 1:])
    elif species == 'nde':
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[:idx])
    elif species == 'aca':
        new_phylogeny.extend([phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[:idx - 1])
    elif species == 'din':
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 2:idx])
        new_phylogeny.extend(phylogeny[:idx - 2])
    elif species == 'fnu':
        new_phylogeny.extend([phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 3:idx - 1])
        new_phylogeny.extend(phylogeny[:idx - 3])
    elif species == 'dap':
        new_phylogeny.extend(phylogeny[idx - 2:idx])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 4:idx - 2])
        new_phylogeny.extend(phylogeny[:idx - 4])
    elif species == 'tid':
        new_phylogeny.extend([phylogeny[idx + 1]])
        new_phylogeny.extend(phylogeny[idx - 3:idx])
        new_phylogeny.extend(phylogeny[idx + 2:])
        new_phylogeny.extend(phylogeny[idx - 5:idx - 3])
        new_phylogeny.extend(phylogeny[:idx - 5])
    elif species == 'aae':
        new_phylogeny.extend([phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx - 4:idx - 1])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 6:idx - 4])
        new_phylogeny.extend(phylogeny[:idx - 6])
    elif species == 'bsu':
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 5:idx])
        new_phylogeny.extend(phylogeny[idx - 7:idx - 5])
        new_phylogeny.extend(phylogeny[:idx - 7])
    elif species == 'lmo':
        new_phylogeny.extend([phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 6:idx - 1])
        new_phylogeny.extend(phylogeny[idx - 8:idx - 6])
        new_phylogeny.extend(phylogeny[:idx - 8])
    elif species == 'sau':
        new_phylogeny.extend(phylogeny[idx - 2:idx])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 7:idx - 2])
        new_phylogeny.extend(phylogeny[idx - 9:idx - 7])
        new_phylogeny.extend(phylogeny[:idx - 9])
    elif species == 'lac':
        new_phylogeny.extend([phylogeny[idx + 1]])
        new_phylogeny.extend(phylogeny[idx - 3:idx])
        new_phylogeny.extend(phylogeny[idx + 2:])
        new_phylogeny.extend(phylogeny[idx - 8:idx - 3])
        new_phylogeny.extend(phylogeny[idx - 10:idx - 8])
        new_phylogeny.extend(phylogeny[:idx - 10])
    elif species == 'snd':
        new_phylogeny.extend([phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx - 4:idx - 1])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 9:idx - 4])
        new_phylogeny.extend(phylogeny[idx - 11:idx - 9])
        new_phylogeny.extend(phylogeny[:idx - 11])
    elif species == 'cpe':
        new_phylogeny.extend(phylogeny[idx - 5:idx])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 10:idx - 5])
        new_phylogeny.extend(phylogeny[idx - 12:idx - 10])
        new_phylogeny.extend(phylogeny[:idx - 12])
    elif species == 'mpn':
        new_phylogeny.extend(phylogeny[idx - 6:idx])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 11:idx - 6])
        new_phylogeny.extend(phylogeny[idx - 13:idx - 11])
        new_phylogeny.extend(phylogeny[:idx - 13])
    elif species == 'syn':
        new_phylogeny.extend(phylogeny[idx + 1:idx + 10])
        new_phylogeny.extend(phylogeny[idx - 7:idx])
        new_phylogeny.extend(phylogeny[idx + 10:])
        new_phylogeny.extend(phylogeny[idx - 12:idx - 7])
        new_phylogeny.extend(phylogeny[idx - 14:idx - 12])
        new_phylogeny.extend(phylogeny[:idx - 14])
    elif species == 'pma':
        new_phylogeny.extend([phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx + 1:idx + 9])
        new_phylogeny.extend(phylogeny[idx - 8:idx - 1])
        new_phylogeny.extend(phylogeny[idx + 9:])
        new_phylogeny.extend(phylogeny[idx - 13:idx - 8])
        new_phylogeny.extend(phylogeny[idx - 15:idx - 13])
        new_phylogeny.extend(phylogeny[:idx - 15])
    elif species == 'cau':
        new_phylogeny.extend(phylogeny[idx - 2:idx])
        new_phylogeny.extend(phylogeny[idx + 1:idx + 8])
        new_phylogeny.extend(phylogeny[idx - 9:idx - 2])
        new_phylogeny.extend(phylogeny[idx + 8:])
        new_phylogeny.extend(phylogeny[idx - 14:idx - 9])
        new_phylogeny.extend(phylogeny[idx - 16:idx - 14])
        new_phylogeny.extend(phylogeny[:idx - 16])
    elif species == 'bbv':
        new_phylogeny.extend(phylogeny[idx + 1:idx + 7])
        new_phylogeny.extend(phylogeny[idx - 3:idx])
        new_phylogeny.extend(phylogeny[idx - 10:idx - 3])
        new_phylogeny.extend(phylogeny[idx + 7:])
        new_phylogeny.extend(phylogeny[idx - 15:idx - 10])
        new_phylogeny.extend(phylogeny[idx - 17:idx - 15])
        new_phylogeny.extend(phylogeny[:idx - 17])
    elif species == 'cgl':
        new_phylogeny.extend([phylogeny[idx + 1], phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx + 2:idx + 6])
        new_phylogeny.extend(phylogeny[idx - 4:idx - 1])
        new_phylogeny.extend(phylogeny[idx - 11:idx - 4])
        new_phylogeny.extend(phylogeny[idx + 6:])
        new_phylogeny.extend(phylogeny[idx - 16:idx - 11])
        new_phylogeny.extend(phylogeny[idx - 18:idx - 16])
        new_phylogeny.extend(phylogeny[:idx - 18])
    elif species == 'mtv':
        new_phylogeny.extend([phylogeny[idx - 1], phylogeny[idx - 2]])
        new_phylogeny.extend(phylogeny[idx + 1:idx + 5])
        new_phylogeny.extend(phylogeny[idx - 5:idx - 2])
        new_phylogeny.extend(phylogeny[idx - 12:idx - 5])
        new_phylogeny.extend(phylogeny[idx + 5:])
        new_phylogeny.extend(phylogeny[idx - 17:idx - 12])
        new_phylogeny.extend(phylogeny[idx - 19:idx - 17])
        new_phylogeny.extend(phylogeny[:idx - 19])
    elif species == 'sco':
        new_phylogeny.extend(phylogeny[idx - 3:idx])
        new_phylogeny.extend(phylogeny[idx + 1:idx + 4])
        new_phylogeny.extend(phylogeny[idx - 6:idx - 3])
        new_phylogeny.extend(phylogeny[idx - 13:idx - 6])
        new_phylogeny.extend(phylogeny[idx + 4:])
        new_phylogeny.extend(phylogeny[idx - 18:idx - 13])
        new_phylogeny.extend(phylogeny[idx - 20:idx - 18])
        new_phylogeny.extend(phylogeny[:idx - 20])
    elif species == 'dra':
        new_phylogeny.extend(phylogeny[idx + 1:idx + 3])
        new_phylogeny.extend(phylogeny[idx - 4:idx])
        new_phylogeny.extend(phylogeny[idx - 7:idx - 4])
        new_phylogeny.extend(phylogeny[idx - 14:idx - 7])
        new_phylogeny.extend(phylogeny[idx + 3:])
        new_phylogeny.extend(phylogeny[idx - 19:idx - 14])
        new_phylogeny.extend(phylogeny[idx - 21:idx - 19])
        new_phylogeny.extend(phylogeny[:idx - 21])
    elif species == 'tth':
        new_phylogeny.extend([phylogeny[idx - 1], phylogeny[idx + 1]])
        new_phylogeny.extend(phylogeny[idx - 5:idx - 1])
        new_phylogeny.extend(phylogeny[idx - 8:idx - 5])
        new_phylogeny.extend(phylogeny[idx - 15:idx - 8])
        new_phylogeny.extend(phylogeny[idx + 2:])
        new_phylogeny.extend(phylogeny[idx - 20:idx - 15])
        new_phylogeny.extend(phylogeny[idx - 22:idx - 20])
        new_phylogeny.extend(phylogeny[:idx - 22])
    elif species == 'fgi':
        new_phylogeny.extend(phylogeny[idx - 2:idx])
        new_phylogeny.extend(phylogeny[idx - 6:idx - 2])
        new_phylogeny.extend(phylogeny[idx - 9:idx - 6])
        new_phylogeny.extend(phylogeny[idx - 16:idx - 9])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 21:idx - 16])
        new_phylogeny.extend(phylogeny[idx - 23:idx - 21])
        new_phylogeny.extend(phylogeny[:idx - 23])
    elif species == 'amo':
        new_phylogeny.extend(phylogeny[idx + 1:idx + 4])
        new_phylogeny.extend(phylogeny[idx - 17:idx])
        new_phylogeny.extend(phylogeny[idx + 4:])
        new_phylogeny.extend(phylogeny[idx - 22:idx - 17])
        new_phylogeny.extend(phylogeny[idx - 24:idx - 22])
        new_phylogeny.extend(phylogeny[:idx - 24])
    elif species == 'tmm':
        new_phylogeny.extend(phylogeny[idx + 1:idx + 3])
        new_phylogeny.extend([phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx - 18:idx - 1])
        new_phylogeny.extend(phylogeny[idx + 3:])
        new_phylogeny.extend(phylogeny[idx - 23:idx - 18])
        new_phylogeny.extend(phylogeny[idx - 25:idx - 23])
        new_phylogeny.extend(phylogeny[:idx - 25])
    elif species == 'cex':
        new_phylogeny.extend([phylogeny[idx + 1]])
        new_phylogeny.extend([phylogeny[idx - 1], phylogeny[idx - 2]])
        new_phylogeny.extend(phylogeny[idx - 19:idx - 2])
        new_phylogeny.extend(phylogeny[idx + 2:])
        new_phylogeny.extend(phylogeny[idx - 24:idx - 19])
        new_phylogeny.extend(phylogeny[idx - 26:idx - 24])
        new_phylogeny.extend(phylogeny[:idx - 26])
    elif species == 'dth':
        new_phylogeny.extend(phylogeny[idx - 3:idx])
        new_phylogeny.reverse()
        new_phylogeny.extend(phylogeny[idx - 20:idx - 3])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 25:idx - 20])
        new_phylogeny.extend(phylogeny[idx - 27:idx - 25])
        new_phylogeny.extend(phylogeny[:idx - 27])
    elif species == 'fsu':
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 21:idx])
        new_phylogeny.extend(phylogeny[idx - 26:idx - 21])
        new_phylogeny.extend(phylogeny[idx - 28:idx - 26])
        new_phylogeny.extend(phylogeny[:idx - 28])
    elif species == 'gau':
        new_phylogeny.extend([phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 22:idx - 1])
        new_phylogeny.extend(phylogeny[idx - 27:idx - 22])
        new_phylogeny.extend(phylogeny[idx - 29:idx - 27])
        new_phylogeny.extend(phylogeny[:idx - 29])
    elif species == 'cph':
        new_phylogeny.extend([phylogeny[idx + 1]])
        new_phylogeny.extend([phylogeny[idx - 2], phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx + 2:])
        new_phylogeny.extend(phylogeny[idx - 23:idx - 2])
        new_phylogeny.extend(phylogeny[idx - 28:idx - 23])
        new_phylogeny.extend(phylogeny[idx - 30:idx - 28])
        new_phylogeny.extend(phylogeny[:idx - 30])
    elif species == 'bfr':
        new_phylogeny.extend([phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx - 3:idx - 1])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 24:idx - 3])
        new_phylogeny.extend(phylogeny[idx - 29:idx - 24])
        new_phylogeny.extend(phylogeny[idx - 31:idx - 29])
        new_phylogeny.extend(phylogeny[:idx - 31])
    elif species == 'rba':
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 4:idx])
        new_phylogeny.extend(phylogeny[idx - 25:idx - 4])
        new_phylogeny.extend(phylogeny[idx - 30:idx - 25])
        new_phylogeny.extend(phylogeny[idx - 32:idx - 30])
        new_phylogeny.extend(phylogeny[:idx - 32])
    elif species == 'cpn':
        new_phylogeny.extend([phylogeny[idx + 1], phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx + 2:])
        new_phylogeny.extend(phylogeny[idx - 5:idx - 1])
        new_phylogeny.extend(phylogeny[idx - 26:idx - 5])
        new_phylogeny.extend(phylogeny[idx - 31:idx - 26])
        new_phylogeny.extend(phylogeny[idx - 33:idx - 31])
        new_phylogeny.extend(phylogeny[:idx - 33])
    elif species == 'ote':
        new_phylogeny.extend([phylogeny[idx - 1], phylogeny[idx - 2]])
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 6:idx - 2])
        new_phylogeny.extend(phylogeny[idx - 27:idx - 6])
        new_phylogeny.extend(phylogeny[idx - 32:idx - 27])
        new_phylogeny.extend(phylogeny[idx - 34:idx - 32])
        new_phylogeny.extend(phylogeny[:idx - 34])
    elif species == 'bbn':
        new_phylogeny.extend(phylogeny[idx + 1:])
        new_phylogeny.extend(phylogeny[idx - 3:idx])
        new_phylogeny.extend(phylogeny[idx - 7:idx - 3])
        new_phylogeny.extend(phylogeny[idx - 28:idx - 7])
        new_phylogeny.extend(phylogeny[idx - 33:idx - 28])
        new_phylogeny.extend(phylogeny[idx - 35:idx - 33])
        new_phylogeny.extend(phylogeny[:idx - 35])
    elif species == 'emi':
        new_phylogeny.extend([phylogeny[idx - 1], phylogeny[idx + 1]])
        new_phylogeny.extend(phylogeny[idx - 4:idx - 1])
        new_phylogeny.extend(phylogeny[idx - 8:idx - 4])
        new_phylogeny.extend(phylogeny[idx - 29:idx - 8])
        new_phylogeny.extend(phylogeny[idx - 34:idx - 29])
        new_phylogeny.extend(phylogeny[idx - 36:idx - 34])
        new_phylogeny.extend(phylogeny[:idx - 36])
    elif species == 'heo':
        new_phylogeny.extend([phylogeny[idx - 2], phylogeny[idx - 1]])
        new_phylogeny.extend(phylogeny[idx - 5:idx - 2])
        new_phylogeny.extend(phylogeny[idx - 9:idx - 5])
        new_phylogeny.extend(phylogeny[idx - 30:idx - 9])
        new_phylogeny.extend(phylogeny[idx - 35:idx - 30])
        new_phylogeny.extend(phylogeny[idx - 37:idx - 35])
        new_phylogeny.extend(phylogeny[:idx - 37])

    return new_phylogeny
