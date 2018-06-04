"""
General interest utilities for the CoMetGeNe pipeline.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import cPickle
import os
import re
import sys


def split_entities(entity_list):
    """Splits a list of entities (such as R numbers or genes) that can contain
    groups of entities (delimited by curly braces) into individual entities.

    :param entity_list: list of entities, possible with groups delimited by
        curly braces
    :return: set of individual entities in entity_list
    """
    separated_entities = set()

    for entity_candidate in entity_list:
        if '{' in entity_candidate:
            entities = re.sub('{', '', entity_candidate)
            entities = re.sub('}', '', entities).split(', ')
            separated_entities = separated_entities.union(entities)
        else:
            separated_entities.add(entity_candidate)

    return separated_entities


def get_abs_dir(directory):
    """Returns the absolute path to the specified directory.

    :param directory: path to a directory
    :return: absolute path to directory
    """
    if not os.path.isdir(directory):
        sys.stderr.write("%s is not a directory.\n" % directory)
        exit(1)
    return os.path.abspath(directory)


def create_directory(directory):
    """Creates the specified directory if it does not exist.

    :param directory: directory to create if it does not exist
    """
    if len(directory) > 0 and not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except OSError:
            raise


def open_device(output):
    """Opens and returns the output device (stdout or a file handler) depending
    on the specified output file. Directories are created as necessary.

    :param output: output file or None if not specified
    :return: stdout if output is not specified or handle to file
    """
    if output is None:
        return sys.stdout
    try:
        create_directory(os.path.dirname(output))
        return open(output, 'w')
    except IOError:
        sys.stderr.write('Cannot open output file \'%s\'.\n' % output)


def pickle(filename, data):
    """Pickles the data to the specified file.

    :param filename: file to pickle the data to
    :param data: data to pickle
    """
    f_pickle = open(filename, 'wb')
    cPickle.dump(data, f_pickle, -1)
    f_pickle.close()


def unpickle(filename):
    """Un-pickles and returns the data from the specified file.

    :param filename: file to un-pickle the data from
    :return: the data contained within the specified file
    """
    f_pickle = open(filename, 'rb')
    data = cPickle.load(f_pickle)
    f_pickle.close()
    return data
