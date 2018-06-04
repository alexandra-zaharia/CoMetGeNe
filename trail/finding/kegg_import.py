"""
Handles data retrieval from KEGG.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import urllib2
import sys
import os
import re
from CoMetGeNeError import CoMetGeNeError
from os.path import exists, isdir, basename

from ..definitions import PICKLE_EC
from ..utils import pickle, unpickle


def check_directory(directory):
    """If the specified directory exists, check whether the permissions are
    sufficient. If the directory does not exist, attempt to create it.

    :param directory: directory where KEGG pathways will be stored in KGML
        format
    """
    if exists(directory):
        if not isdir(directory):
            raise CoMetGeNeError(
                'import_not_dir', basename(__file__), directory)
        else:
            if not os.access(directory, os.R_OK):
                raise CoMetGeNeError(
                    'import_not_r', basename(__file__), directory)
            if not os.access(directory, os.W_OK):
                raise CoMetGeNeError(
                    'import_not_w', basename(__file__), directory)
            if not os.access(directory, os.X_OK):
                raise CoMetGeNeError(
                    'import_not_x', basename(__file__), directory)
    else:  # Create the output directory if it doesn't exist.
        try:
            os.makedirs(directory)
        except Exception:
            raise CoMetGeNeError('import_mkdir', basename(__file__), directory)


def remove_meta_pathways(pathways):
    """Removes meta-pathway maps, i.e. maps with an ID starting from 01100, from
    a list of KEGG pathway map IDs.

    The 'pathways' list is modified.

    :param pathways: list of strings designating KEGG pathway map IDs
    """
    ids_to_remove = list()
    for pathway in pathways:
        digit_match = re.search('\d{5}', pathway)
        assert digit_match
        kegg_id = int(digit_match.group())
        if kegg_id >= 1100:
            ids_to_remove.append(pathway)
    for meta_pw in ids_to_remove:
        pathways.remove(meta_pw)


def save_pathways(pathways, organism, directory):
    """Downloads the specified pathway maps from KEGG in KGML format for a given
    species and saves them to the specified directory.

    If the output file already exists, it is not overwritten.

    :param pathways: list of KEGG pathway map IDs
    :param organism: species for which to retrieve the pathways
    :param directory: directory where pathway maps for 'organism' will be stored
        in KGML format
    """
    saved_files = 0
    for i in range(len(pathways)):
        kgml = os.path.join(directory, 'path_' + pathways[i] + '.kgml')

        if exists(kgml):
            print "\tOutput file %s already exists. Skipping." % kgml
            sys.stdout.flush()
        else:
            print "\tRetrieving pathway", pathways[i], \
                  "(" + str(i+1) + "/" + str(len(pathways)) + ") ...",
            sys.stdout.flush()
            pathway_url = 'http://rest.kegg.jp/get/' + pathways[i] + '/kgml'
            try:
                data = urllib2.urlopen(pathway_url)
                xml = data.read()
                print "done"

                print "\t\tWriting pathway to file", kgml, "...",
                sys.stdout.flush()
                k_out = open(kgml, 'w')
                k_out.write(xml)
                k_out.close()
                saved_files += 1
                print "done"
            except urllib2.URLError:
                fmt_msg = "\n%s: HTTP download error\n" % \
                          os.path.basename(__file__)
                sys.stderr.write(fmt_msg)
            except IOError:
                fmt_msg = "\n%s: Error writing to file\n" % \
                          os.path.basename(__file__)
                sys.stderr.write(fmt_msg)
    print "Done! Saved %d pathway(s) for '%s' under %s.\n" % (
        saved_files, organism, directory)


def download_kgml(organism, directory):
    """Downloads all non-meta-pathways for a given species and saves them to the
    specified directory unless the output files already exist.

    A KEGG meta-pathway has an ID starting from 01100.

    :param organism: species for which to retrieve the pathways
    :param directory: directory where pathway maps for 'organism' will be stored
        in KGML format
    """
    if len(organism) != 3 and len(organism) != 4:
        raise CoMetGeNeError('import_org_code', basename(__file__), organism)

    # Ensure there are sufficient permissions on the output directory, and
    # create it if necessary.
    check_directory(directory)

    # Retrieve a list of metabolic pathways for the query organism.
    print "Retrieving metabolic pathways for '" + organism + "':"
    sys.stdout.flush()
    pathway_list_url = 'http://rest.kegg.jp/list/pathway/' + organism
    try:
        data = urllib2.urlopen(pathway_list_url)
        pathways = data.read().split('\n')
        pathways = filter(lambda x: x != '', pathways)  # remove empty lines
        pathways = [pw.split()[0].split(':')[1] for pw in pathways]
        remove_meta_pathways(pathways)
    except urllib2.URLError:
        raise CoMetGeNeError('import_not_found', basename(__file__), organism)

    if len(pathways) == 0:
        raise CoMetGeNeError('import_not_found', basename(__file__), organism)

    save_pathways(pathways, organism, directory)


def retrieve_ec_numbers():
    """Retrieves and returns EC number and R number associations, downloading
    the required information from KEGG if necessary.

    If the destination file designated by PICKLE_EC already exists, its contents
    is simply loaded (un-pickled) and returned.

    :return: dict associating R numbers (keys) to EC numbers (values)
    """
    if not os.path.exists(PICKLE_EC):
        query_url = 'http://rest.kegg.jp/link/reaction/enzyme'
        print "Retrieving EC numbers from KEGG ...",
        try:
            data = urllib2.urlopen(query_url)
            results = data.read().split('\n')
        except urllib2.URLError:
            raise CoMetGeNeError('import_ec', basename(__file__), None)

        ec_numbers = dict()
        for line in results:
            if len(line) > 0:  # Ignore empty lines.
                ec_data = line.split('\t')
                assert len(ec_data) == 2
                ec = ec_data[0].replace('ec:', '')
                reaction = ec_data[1]
                if reaction not in ec_numbers:
                    ec_numbers[reaction] = list()
                ec_numbers[reaction].append(ec)
        print "done\n"

        pickle(PICKLE_EC, ec_numbers)
    else:
        ec_numbers = unpickle(PICKLE_EC)

    return ec_numbers


def get_slices(data, slice_size):
    """Slices up and returns the data in slices of slice_size.

    :param data: list to divide in one or several slices of size slice_size
    :param slice_size: integer designating the size of a slice from data
    :return: list of len(data) / slice_size slices of data of size slice_size if
        the number of items in data is a multiple of slice_size, or list of
        len(data) / slice_size + 1 slices of data of size slice_size except for
        the last slice, of size len(data) - slice_size * len(data) / slice_size
    """
    slices = list()

    indexes = [i for i in range(0, len(data), slice_size)]

    for i in range(0, len(indexes) - 1):
        slices.append(data[indexes[i]:indexes[i + 1]])

    if len(data) > indexes[-1]:  # is there a last slice?
        slices.append(data[indexes[-1]:])

    return slices


def extract_gene_info(gene_info, org, genomes):
    """Stores information on a given protein-coding gene of species 'org' in the
    dict 'genomes'.

    The 'genomes' dict is modified unless 'gene_info' does not designate a CDS.

    :param gene_info: textual information on a gene entry for species 'org', as
        retrieved from KEGG GENES
    :param org: species to which belongs the gene description gene_info
    :param genomes: dict of dicts storing gene information for every gene of
        every species in the dict, namely the name of the chromosome on which
        the gene is located, the strand on the chromosome, as well as the
        position of the gene on the chromosome (in nucleotides)
    """
    if gene_info[0].split()[2] != 'CDS':
        return

    gene = org + ':' + gene_info[0].split()[1]
    genomes[org][gene] = dict()

    pos_info = ''
    for entry in gene_info:
        if entry.split()[0] == 'POSITION':
            pos_info = entry.split()[1]
            break
    assert len(pos_info) > 0

    fields = pos_info.split(':')
    if len(fields) == 1:  # Only one chromosome.
        genomes[org][gene]['chr'] = 'chromosome'
        pos = 0
    else:
        genomes[org][gene]['chr'] = 'chromosome ' + fields[0]
        pos = 1
    genomes[org][gene]['fwd'] = False if 'complement' in fields[pos] else True

    genomes[org][gene]['pos'] = list()
    positions = fields[pos].split(',')
    for position in positions:
        pos_match = re.search('(\d+).*\.\..*(\d+)', position)
        assert pos_match
        start, end = pos_match.group().split('..')
        start = int(re.search('(\d+)', start).group())
        end = int(re.search('(\d+)', end).group())
        genomes[org][gene]['pos'].append((start, end))


def retrieve_gene_info(genes, organism, genomes):
    """For every KEGG GENES entry in 'genes' designating multiple genes of
    species 'organism', retrieves the relevant gene information for coding
    sequences and stores it in the dict 'genomes'.

    The 'genomes' dict is modified if at least one gene in 'genes' is a CDS.

    :param genes: list of strings designating multiple KEGG GENES entries
        separated by a line containing '///'
    :param organism: species to which the gene entries in 'genes' belong to
    :param genomes: dict of dicts storing gene information for every gene of
        every species in the dict, namely the name of the chromosome on which
        the gene is located, the strand on the chromosome, as well as the
        position of the gene on the chromosome (in nucleotides)
    """
    query_url = 'http://rest.kegg.jp/get/' + '+'.join(genes)
    data = urllib2.urlopen(query_url).read().split('\n')

    indices = [i for i, x in enumerate(data) if x == "///"]

    if len(indices) == 0:
        extract_gene_info(data, organism, genomes)
    else:
        for i in range(len(indices) - 1):
            if i == 0:
                extract_gene_info(data[:indices[i]], organism, genomes)
            extract_gene_info(
                data[(indices[i] + 1):indices[i + 1]], organism, genomes)


def retrieve_genome_info(organism, genomes=None):
    """If species 'organism' is not present in the dict 'genomes' storing gene
    information for several species, retrieves all genomic information from
    KEGG GENES and stores it in the dict 'genomes'.

    If 'genomes' is not None and 'organism' is not a key of 'genomes', then the
    dict 'genomes' is modified.

    The return value of this method is retrieved by CoMetGeNe.py, and not
    retrieved by CoMetGeNe_launcher.py, respectively.

    :param organism: species for which genomic information will be retrieved if
        not already present in the dict 'genomes'
    :param genomes: dict of dicts storing gene information for every gene of
        every species in the dict, namely the name of the chromosome on which
        the gene is located, the strand on the chromosome, as well as the
        position of the gene on the chromosome (in nucleotides)
    :return: the species for which genomic information has been retrieved
        ('organism') and the associated genomic information
    """
    if genomes is not None and organism in genomes:
        return

    genes_dict = genomes if genomes is not None else dict()
    genes_dict[organism] = dict()
    list_url = 'http://rest.kegg.jp/list/' + organism
    kegg_gene_ids = list()

    try:
        data = urllib2.urlopen(list_url).read().split('\n')
        data = filter(lambda x: x != '', data)  # Remove empty lines.
        for line in data:
            kegg_gene_ids.append(line.split()[0])

        slice_size = 10
        gene_counter = 0
        for genes in get_slices(kegg_gene_ids, slice_size):
            gene_counter = min(gene_counter + slice_size, len(kegg_gene_ids))
            progress = "Retrieving information for '%s' genes: %d/%d" % (
                organism, gene_counter, len(kegg_gene_ids))
            sys.stdout.write('%s\r' % progress)
            sys.stdout.flush()
            retrieve_gene_info(genes, organism, genes_dict)

        sys.stdout.write('\n')
    except urllib2.URLError:
        fmt_msg = "\n%s: HTTP download error\n" % os.path.basename(__file__)
        sys.stderr.write(fmt_msg)

    return organism, genes_dict[organism]
