#!/usr/bin/env python2

"""
This script performs trail finding for a given species.

Retrieval of metabolic pathways, genomic information and EC numbers associations
is also handled, as needed.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import multiprocessing
import os
import sys
import time

import trail.finding.consistency as consistency
import trail.finding.graph as graph
import trail.finding.HNet as HNet
import trail.finding.kegg_import as kegg_import
from trail.finding.Exclusions import Exclusions
from trail.finding.CoMetGeNeError import CoMetGeNeError, error
from trail.finding.NetworkBuilder import NetworkBuilder
from trail.finding.output import output_trail
from trail.utils import open_device
from parsers import arg_parser, kgml_parser


def run_HNet(kgml, args, exclusions, network_inst, dev_out):
    """Runs HNet (trail finding) for a given pathway of a given species.

    If trail finding takes too long, the analysis is aborted and the pathway is
    blacklisted for the gap parameters with which CoMetGeNe was executed.

    :param kgml: metabolic pathway in KGML format
    :param args: command-line arguments for this script
    :param exclusions: Exclusions object representing blacklisted pathways
    :param network_inst: NetworkBuilder object
    :param dev_out: output device for results (file or stdout)
    :return: list of trails found by the HNet algorithm for the given pathway
    """
    # The HNet_on_every_arc function of HNet.py is started as a process that is
    # terminated if it takes longer than a set timeout. Terminating the process
    # also results in blacklisting the pathway, as its analysis takes too long.
    HNet_queue = multiprocessing.Queue()
    HNet_process = multiprocessing.Process(
        target=HNet.HNet_on_every_arc,
        args=(HNet_queue, network_inst.G, network_inst.D,
              network_inst.reactions,)
    )
    HNet_process.start()
    HNet_process.join(timeout=args.timeout)
    if HNet_process.is_alive():
        HNet_process.terminate()
        HNet_process.join()
        aborted = kgml + ': Aborted (analysis takes longer than '
        aborted += str(args.timeout) + ' seconds)\n\n'
        dev_out.write(aborted)
        exclusions.blacklist(
            network_inst, 
            os.path.join(os.path.abspath(args.DIR), kgml))
        return None

    return HNet_queue.get_nowait()


def check_consistency(trail, network_inst):
    """Checks whether the given trail is consistent with the original reaction
    and gene networks.

    The trail is a solution for a metabolic pathway and an undirected graph
    built on the same vertex set as the metabolic pathway, representing gene
    neighborhood in terms of reactions. This method tests whether this trail is
    still a solution for the metabolic pathway and the original undirected graph
    representing gene neighborhood.

    If the consistency check takes longer than 30 seconds, it is aborted.

    :param trail: a trail produced by the HNet algorithm
    :param network_inst: NetworkBuilder object
    :return: True if the trail is consistent with the original reaction and gene
        networks, False if it is not, and None if the consistency check takes
        longer than 30 seconds
    """
    cons_queue = multiprocessing.Queue()
    cons_process = multiprocessing.Process(
        target=consistency.is_consistent,
        args=(cons_queue, network_inst.G_reduced, network_inst.reactions,
              trail,)
    )
    cons_process.start()
    cons_process.join(timeout=30)
    if cons_process.is_alive():
        cons_process.terminate()
        cons_process.join()
        return None

    return cons_queue.get_nowait()


def analyze_kgml(kgml, args, G_init, exclusions, ec_numbers, dev_out):
    """Runs the HNet algorithm on the given KGML file (metabolic pathway) and
    outputs results.

    :param kgml: metabolic pathway in KGML format
    :param args: command-line arguments for this script
    :param G_init: undirected graph representing gene neighborhood for the
        given species
    :param exclusions: Exclusions object representing blacklisted pathways
    :param ec_numbers: dict associating a list of EC numbers (values) to R
        numbers (keys)
    :param dev_out: output device for results (file or stdout)
    """
    pathway = os.path.join(os.path.abspath(args.DIR), kgml)
    if not exclusions.can_analyze(pathway, args.delta_G, args.delta_D):
        return

    title = kgml + ': ' + kgml_parser.get_pathway_title(pathway)
    dev_out.write(title + '\n')

    network_inst = NetworkBuilder(
        G_init, kgml, args, exclusions, ec_numbers, dev_out)

    # Proceed only if all networks have been initialized within the allotted
    # timeout.
    if not network_inst.blacklisted:
        trails = run_HNet(kgml, args, exclusions, network_inst, dev_out)
        if trails is not None:
            if len(trails) == 0:
                dev_out.write(kgml + ': (not found)\n')
            else:
                for trail in trails:
                    if check_consistency(trail, network_inst):
                        output_trail(
                            kgml,
                            trail,
                            network_inst,
                            dev_out)
                dev_out.write(kgml + ':\n')

            elapsed = '%s: --- %.2f seconds ---\n\n' % (
                kgml, (time.time() - network_inst.start_pw))
            dev_out.write(elapsed)


def main():
    """Performs trail finding (HNet) for a given species.

    Metabolic pathways for the given species are retrieved if necessary, as well
    as its genomic information.

    Results are either stored in a file, or displayed on stdout.
    """
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)  # unbuffered mode
    start_time = time.time()  # record the starting time
    args = arg_parser.parse_cmd_arguments()
    dev_out = open_device(args.output)

    try:
        if not args.skip_import:
            kegg_import.download_kgml(args.ORG, args.DIR)

        deltas = '--- delta_G = %d, delta_D = %d ---\n\n' % (
            args.delta_G, args.delta_D)
        dev_out.write(deltas)

        # Build undirected graph representing the genome.
        G_init = graph.build_undirected_graph(args)

        # Determine which pathways should be skipped for the current species.
        exclusions = Exclusions(args.ORG)

        # Retrieve the associations between EC numbers and reactions.
        ec_numbers = kegg_import.retrieve_ec_numbers()

        # Determine the list of kgml files to analyze.
        directory = os.path.abspath(args.DIR)
        pathways = [filename for filename in os.listdir(directory)
                    if filename.lower().endswith('.kgml')]
        pathways.sort()
        
        # Run CoMetGeNe for every metabolic pathway in the specified directory.
        for kgml in pathways:
            analyze_kgml(kgml, args, G_init, exclusions, ec_numbers, dev_out)

        elapsed = '--- %.2f seconds ---\n' % (time.time() - start_time)
        dev_out.write(elapsed)

        if args.output is not None:
            dev_out.close()

    except CoMetGeNeError as err:
        sys.stderr.write(err.text + '\n')
        if args.output is not None:
            dev_out.close()
        exit(error[err.value])
    else:
        exit(0)


if __name__ == '__main__':
    main()
