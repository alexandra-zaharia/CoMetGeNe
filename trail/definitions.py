"""
These definitions are needed for both trail finding and trail grouping.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import os

EXCLUSIONS_FILENAME = 'excluded_pathways.txt'

PICKLE_DIR = 'pickle'  # directory storing .pickle files

PICKLE_EC = os.path.join(PICKLE_DIR, 'kegg_ec_info.pickle')

PICKLE_RN_FILENAME = 'grouping_r_numbers.pickle'
PICKLE_RN = os.path.join(PICKLE_DIR, PICKLE_RN_FILENAME)

PICKLE_RS_FILENAME = 'grouping_reaction_sets.pickle'
PICKLE_RS = os.path.join(PICKLE_DIR, PICKLE_RS_FILENAME)

PICKLE_GEN_FILENAME = 'grouping_species_object.pickle'
PICKLE_GEN = os.path.join(PICKLE_DIR, PICKLE_GEN_FILENAME)

PICKLE_GENOME_FILENAME = 'kegg_genome_info.pickle'
PICKLE_GENOME = os.path.join(PICKLE_DIR, PICKLE_GENOME_FILENAME)

MAX_DELTA = 3  # maximum number of skipped genes
