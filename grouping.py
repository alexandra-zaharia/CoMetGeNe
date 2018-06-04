#!/usr/bin/env python2

"""
This script performs trail grouping for a given reference species.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""
import argparse

import sys

import os

from trail.grouping.TrailGroupingError import TrailGroupingError
from parsers.CoMetGeNe_parser import parse_CoMetGeNe
from trail.grouping.phylogeny import warning_phylogeny
from trail.grouping.trail_grouping_by_genes import fill_gene_matrix, \
    print_species_table
from trail.grouping.trail_grouping_by_reactions import get_symbols, \
    print_symbols
from trail.grouping.trail_grouping_utils import initialization
from trail.utils import open_device


def parse_cmd_arguments():
    """Parses command-line arguments for trail grouping.

    :return: command-line arguments for trail grouping
    """
    desc = 'Groups CoMetGeNe trails by either genes or reactions, producing ' \
        'a CSV file.'
    example = \
        '''
KGML needs to contain a subdirectory for every species for which a result file 
is present in RESULTS. The subdirectory names need to be the three- or four-
letter KEGG codes for the species in question (e.g., \'bsu\', \'eco\', \'pae\', 
etc.). Each species subdirectory is expected to contain metabolic pathways in
KGML format.  

Example: running 

    python2 grouping.py genes results/ data/ eco -o grouping_gene_eco.csv
    
will perform trail grouping by genes for the reference species 'eco'. The 
CoMetGeNe results are stored in 'results/', and the KGML files are available in
'data/'.
A CSV file is produced ('grouping_gene_eco.csv').
        '''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=desc,
        epilog=example)

    parser._positionals.title = "Required arguments"
    parser._optionals.title = "Optional arguments"

    parser.add_argument(
        'type',
        choices=['genes', 'reactions'],
        help='type of trail grouping to perform '
             '(possible values: \'genes\' or \'reactions\')'
    )

    parser.add_argument(
        'RESULTS', help='directory storing CoMetGeNe results')

    parser.add_argument(
        'KGML', help='directory containing input KGML files'
    )

    parser.add_argument(
        'ORG', help='reference species (KEGG organism code)')

    parser.add_argument(
        '--output', '-o', help='output file (CSV)')

    return parser.parse_args()


def main():
    """Performs trail grouping for a given reference species, by genes or by
    reactions.
    """
    args = parse_cmd_arguments()
    dev_out = open_device(args.output)
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    warning_phylogeny(args.KGML)

    print 'Parsing CoMetGeNe results ...',
    results = parse_CoMetGeNe(args.RESULTS)
    print 'done'
    try:
        if args.type == 'genes':
            species = fill_gene_matrix(results, args.ORG, args.KGML)
            print_species_table(species, args.KGML, dev_out)
        else:
            init_data_structs = initialization(
                results, args.ORG, args.KGML, True)
            symbols = get_symbols(init_data_structs, args.ORG, args.KGML)
            print_symbols(
                init_data_structs, symbols, args.ORG, args.KGML, dev_out)
    except TrailGroupingError as tge:
        tge.display_error()
        sys.exit(1)


if __name__ == "__main__":
    main()
