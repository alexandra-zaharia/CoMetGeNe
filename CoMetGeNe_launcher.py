#!/usr/bin/env python2

"""
This script allows parallel execution of CoMetGeNe.py, especially useful when
performing trail finding for a large number of species.

Version: 1.1 (October 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""
import os
import subprocess
import multiprocessing

from trail.definitions import PICKLE_GENOME
from trail.finding.kegg_import import download_kgml, retrieve_genome_info, \
    retrieve_ec_numbers
from trail.utils import create_directory, unpickle, pickle

################################################################################
#           Configure the following variables to your liking
################################################################################

# Add three- or four-letter KEGG codes to this list. CoMetGeNe will be executed
# on these species.
org_codes = [
    'aae', 'aca', 'afi', 'amo', 'ara', 'bbn', 'bbv', 'bfr', 'bsu', 'cau',
    'cex', 'cgl', 'cpe', 'cph', 'cpn', 'dap', 'din', 'dra', 'dth', 'eco',
    'emi', 'fgi', 'fnu', 'fsu', 'gau', 'gsu', 'heo', 'lac', 'lmo', 'mpn',
    'mtv', 'nde', 'nme', 'ote', 'pae', 'pma', 'rba', 'rrj', 'rso', 'sau',
    'sco', 'snd', 'spc', 'syn', 'tid', 'tmm', 'tth', 'vco', 'xfa', 'ype'
]

delta_genes_max = 3  # maximum number of genes that can be skipped
delta_reactions_max = 3  # maximum number of reactions that can be skipped

# This is the directory where metabolic pathways are stored. For each of the
# species in the list 'org_codes', a subdirectory with the same name (three- or
# four-letter KEGG code) either exists in 'kgml_dir', or is created. Each of
# these subdirectories either contains metabolic pathways for the species in
# question or will be storing them once CoMetGeNe is executed.
kgml_dir = 'data'

results_dir = 'results'  # directory where CoMetGeNe results are stored

# This is the maximum number of simultaneous connections to KEGG for pathway
# retrieval. Decrease this value if pathway retrieval fails. As of April 2018,
# it would appear that 3 simultaneous connections are accepted.
kegg_max_thr_pw = 3

# This is the maximum number of simultaneous connections to KEGG for genomic
# information retrieval. Decrease this value if genomic information retrieval
# fails. As of April 2018, it would appear that 2 simultaneous connections are
# accepted.
kegg_max_thr_gen = 2

# This is the number of threads on which CoMetGeNe is ran (by default, the 
# maximum number of available threads).
n_thr_cometgene = multiprocessing.cpu_count()

################################################################################
#               CoMetGeNe launcher script for parallel execution
################################################################################
def retrieve_pathways():
    """Retrieves metabolic pathways from KEGG for every species in the data set.

    As of May 2018, the maximum number of threads for pathway retrieval is 2.
    """
    pool = multiprocessing.Pool(
        min(kegg_max_thr_pw, len(org_codes), multiprocessing.cpu_count()))

    for org in org_codes:
        org_dir = os.path.join(kgml_dir, org)
        if not os.path.exists(org_dir) or len(os.listdir(org_dir)) == 0:
            pool.apply_async(download_kgml, (org, org_dir,))

    pool.close()
    pool.join()


def retrieve_genomes(genomes):
    """Retrieves genomic information (gene names, chromosome, strand, position)
    from KEGG for every species in the data set.

    As of May 2018, the maximum number of threads for genomic information
    retrieval is 3.

    :param genomes: data structure containing genomic information.
    :return:
    """
    pool = multiprocessing.Pool(
        min(kegg_max_thr_gen, len(org_codes), multiprocessing.cpu_count()))
    manager = multiprocessing.Manager()
    lock = manager.Lock()

    for org in org_codes:
        if org not in genomes:
            pool.apply_async(retrieve_genome_info, (org, None, lock,))

    pool.close()
    pool.join()


def run_CoMetGeNe():
    """Runs CoMetGeNe.py on all available threads for every species in the
    data set and for every combinations of the gap parameters.

    The gap parameters specify how many genes (delta_G) and reactions (delta_D)
    can be skipped.
    """
    pool = multiprocessing.Pool(n_thr_cometgene)
    create_directory(results_dir)

    for delta_G in range(0, delta_genes_max + 1):
        for delta_D in range(0, delta_reactions_max + 1):
            for org in org_codes:
                output = os.path.join(
                    results_dir,
                    org + '_dG' + str(delta_G) + '_dD' + str(delta_D) + '.hnet'
                )
                org_dir = os.path.join(kgml_dir, org)

                CoMetGeNe = [
                    'python2', 'CoMetGeNe.py',
                    '-s', org, org_dir,
                    '-dG', str(delta_G),
                    '-dD', str(delta_D),
                    '-o', output]
                print ' '.join(CoMetGeNe)
                pool.apply_async(subprocess.call, (CoMetGeNe,))

    pool.close()
    pool.join()


def main():
    """Runs CoMetGeNe.py for every species in the data set.

    For every species in the data set, metabolic pathways and genomic
    information are retrieved from KEGG. Then EC numbers associations are also
    extracted. Finally, CoMetGeNe.py is ran for all species, for all
    combinations of the gap parameters.

    Output files and directories are created as needed.
    """
    retrieve_pathways()

    if not os.path.exists(PICKLE_GENOME):
        genomes = dict()
    else:
        genomes = unpickle(PICKLE_GENOME)

    retrieve_genomes(genomes)
    
    retrieve_ec_numbers()

    run_CoMetGeNe()


if __name__ == '__main__':
    main()

