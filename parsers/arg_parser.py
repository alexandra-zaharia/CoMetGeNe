"""
This is the command-line parser for CoMetGeNe.py.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import argparse


def parse_cmd_arguments():
    """Parses command-line arguments for CoMetGeNe.py.

    :return: command-line arguments for CoMetGeNe
    """
    desc = \
        '''
Determines maximum trails of reactions for the specified organisms such that the 
genes encoding the enzymes involved in the trails are neighbors.

A trail of reactions is a sequence of reactions that can repeat reactions 
(vertices), but not arcs between reactions.

Metabolic pathways and genomic information are automatically retrieved from the 
KEGG knowledge base.
        '''
    example = \
        '''
Example: running

    python2 CoMetGeNe.py eco data/ -dG 2 -o eco.out
    
downloads metabolic pathways for species 'eco' to directory 'data/'. Trail 
finding is performed, allowing two genes to be skipped at most (-dG 2). 
Reactions cannot be skipped (-dD is 0 by default). Maximum trails of reactions 
such that the reactions are catalyzed by products of neighboring genes are saved 
in the output file 'eco.out'.
        '''
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=desc,
        epilog=example
    )

    parser._positionals.title = "Required arguments"
    parser._optionals.title = "Optional arguments"

    parser.add_argument(
        'ORG',
        help="query organism (three- or four-letter KEGG code, e.g. 'eco' " +
             "for Escherichia coli K-12 MG1655). See full list of KEGG " +
             "organism codes at http://rest.kegg.jp/list/genome"
    )

    parser.add_argument(
        'DIR',
        help="directory storing metabolic pathways for the query organism " +
             "ORG or where metabolic pathways for ORG will be downloaded"
    )

    parser.add_argument(
        '--delta_G', '-dG',
        type=get_delta_value,
        metavar='NUMBER',
        default=0,
        help="the NUMBER of genes that can be skipped (default: 0)"
    )

    parser.add_argument(
        '--delta_D', '-dD',
        type=get_delta_value,
        metavar='NUMBER',
        default=0,
        help="the NUMBER of reactions that can be skipped (default: 0)"
    )

    parser.add_argument(
        '--timeout', '-t',
        type=get_timeout,
        metavar='SECONDS',
        default=300,
        help="timeout in SECONDS (default: 300)"
    )

    parser.add_argument(
        '--output', '-o',
        help='output file'
    )

    parser.add_argument(
        '--skip-import', '-s',
        action='store_true',
        help="skips importing metabolic pathways from KEGG, attempting to " +
             "use locally stored KGML files if they are present under the " +
             "specified directory (DIR)"
    )

    parser.add_argument(
        '--both-strands', '-b',
        action='store_true',
        help="considers neighboring genes on both strands of a given " +
             "chromosome (by default, only genes located on a single strand " +
             "are considered neighbors)"
    )

    return parser.parse_args()


def get_delta_value(value):
    """Returns the value associated to the delta parameters as an integer, and
    ensures that this value it is between 0 and 10.

    :param value: the value assigned to either of the delta parameters
    :return: the integer value assigned to either of the delta parameters
    """
    if not value.isdigit():
        raise argparse.ArgumentTypeError("contains non-digit characters.")

    delta = int(value)
    if delta < 0 or delta > 10:
        raise argparse.ArgumentTypeError("must be between 0 and 10.")

    return delta


def get_timeout(value):
    """Returns the CoMetGeNe timeout (in seconds).

    :param value: the value assigned to the timeout
    :return: the integer value assigned to the timeout
    """
    if not value.isdigit():
        raise argparse.ArgumentTypeError("contains non-digit characters.")

    timeout = int(value)
    if timeout < 10:
        raise argparse.ArgumentTypeError("must be at least 10.")

    return timeout
