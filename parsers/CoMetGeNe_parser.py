"""
Utilities for parsing CoMetGeNe results.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import os

import re

from Trails import ReactionTrails
from trail.utils import get_abs_dir


def get_organism(filename):
    """Returns the KEGG organism code for the CoMetGeNe result file 'filename'.

    :param filename: CoMetGeNe result file
    :return: KEGG organism code
    """
    organism = None
    pattern = re.compile("path_([a-z]{3,4})")
    with open(filename, 'r') as f_in:
        for line in f_in:
            if line.startswith('path_'):
                org = pattern.search(line)
                assert org
                if organism is None:
                    organism = org.groups()[0]
                else:
                    if organism != org.groups()[0]:
                        return None
    return organism


def get_max_deltas(filename):
    """Returns maximum gap parameters for the specified CoMetGeNe result file.

    :param filename: CoMetGeNe result file
    :return: tuple representing the gap parameters for genome and metabolic
        pathway, respectively
    """
    pattern = re.compile("--- delta_G = (\d+), delta_D = (\d+).* ---")
    with open(filename, 'r') as f_in:
        for line in f_in:
            deltas = pattern.search(line)
            if deltas:
                dG, dD = deltas.groups()
                return int(dG), int(dD)
    return None, None


def parse_entities(line):
    """Parses a line for R numbers, EC numbers, or KGML internal reaction IDs.

    :param line: line in CoMetGeNe result file
    :return: list of reaction identifiers (KEGG R numbers), EC numbers, or KGML
        internal reaction IDs
    """
    entities = line.split(':')[1].strip()
    entities = re.sub("\[.*?\]", "", entities)  # remove skipped entities
    entities = re.sub("-> {2}->", "->", entities)  # remove double arrows
    return entities.split(' -> ')


def get_trails(filename, org, results):
    """Parses CoMetGeNe result file for the specified organism and adds new
    trails accordingly into the data structure used for storing results.

    The ReactionTrails results object is modified.

    :param filename: CoMetGeNe result file
    :param org: KEGG organism code
    :param results: ReactionTrails object storing CoMetGeNe trails and related
        information
    """
    pw_pattern = re.compile('path_' + org + '(\d{5})')
    r_pattern = re.compile('-> R(\d){5}')
    ec_pattern = re.compile('(\d)(\.\d{1,3}){3}|\?')
    rid_pattern = re.compile('path_' + org + '(\d{5})' + '.kgml: ' + '(\d+)')
    max_dG, max_dD = get_max_deltas(filename)
    with open(filename, 'r') as f_in:
        reactions = None
        ec_numbers = None
        genes = None
        r_ids = None
        for line in f_in:
            if r_pattern.search(line):
                reactions = parse_entities(line)
            elif ec_pattern.search(line):
                ec_numbers = parse_entities(line)
            elif rid_pattern.search(line):
                r_ids = parse_entities(line)
            elif re.search(org + ':', line) and '->' in line:
                genes = ' '.join(line.split()[1:])
                # Remove skipped genes.
                genes = re.sub("\[.*?\]", "", genes)
                # Remove double arrows.
                genes = re.sub("-> {2}->", "->", genes)
                genes = genes.split(' -> ')
            else:
                pathway = pw_pattern.match(line)
                if pathway:
                    p_id = pathway.groups()[0]
                    if line == 'path_' + org + p_id + '.kgml:\n' \
                            and genes is not None:
                        if 'ec_numbers' not in vars():
                            ec_numbers = list()
                        trail = results.add_trail(reactions, ec_numbers)
                        assert trail is not None
                        trail.add_organism(org)
                        trail.add_pathway(org, p_id, r_ids, genes)
                        organism = trail.find_by_o_id(org)
                        assert organism is not None
                        pw = organism.find_by_p_id(p_id)
                        assert pw is not None
                        ir = pw.find_by_r_id(r_ids)
                        assert ir is not None
                        ir.set_max_deltas(max_dG, max_dD)


def parse_CoMetGeNe(directory):
    """Parses CoMetGeNe results in the specified directory and returns the
    reaction trails as a ReactionTrails object.

    :param directory: path to the directory containing CoMetGeNe results
    :return: ReactionTrails object storing parsed CoMetGeNe results
    """
    CoMetGeNe_directory = get_abs_dir(directory)

    files = list()
    for CoMetGeNe_output in os.listdir(CoMetGeNe_directory):
        filename = os.path.join(CoMetGeNe_directory, CoMetGeNe_output)
        if os.path.isfile(filename):
            files.append(filename)

    results = ReactionTrails()

    for filename in files:
        org = get_organism(filename)
        assert org is not None, \
            'Organism could not be determined for %s' % filename
        get_trails(filename, org, results)

    return results
