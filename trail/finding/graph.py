"""
These utilities handle graphs.

Version: 1.1 (October 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import os
import networkx as nx

from kegg_import import retrieve_genome_info
from ..definitions import PICKLE_GENOME
from ..utils import unpickle, pickle


def build_undirected_graph(args):
    """Builds and returns the gene neighborhood graph for the specified species.

    If no genomic information is available for the designated species, it will
    be downloaded from KEGG.

    :param args: command-line arguments for CoMetGeNe.py
    :return: undirected graph representing gene neighborhood for the specified
        species
    """
    if not os.path.exists(PICKLE_GENOME):
        genomes = dict()
    else:
        genomes = unpickle(PICKLE_GENOME)

    if args.ORG not in genomes:
        retrieve_genome_info(args.ORG, genomes)
        pickle(PICKLE_GENOME, genomes)

    chr_list = set()
    for gene_id in genomes[args.ORG]:
        chr_list.add(genomes[args.ORG][gene_id]['chr'])
    return build_gene_graph(genomes[args.ORG], list(chr_list),
                            args.both_strands)


def build_gene_graph(gene_data, chromosomes, both_strands=False):
    """Builds and returns an undirected graph representing gene neighborhood.

    :param gene_data: dict storing chromosome name, strand, and gene position,
        for every gene in the dict
    :param chromosomes: list of chromosome names
    :param both_strands: bool specifying whether genes on both strands of a
        given chromosome should be considered neighbors (if set to True) or only
        genes on a single strand should be considered neighbors (if set to
        False)
    :return: undirected graph representing gene neighborhood
    """

    G = nx.Graph()
    for gene in gene_data:
        G.add_node(gene)

    for chromosome in chromosomes:
        genes = _get_genes_on_chromosome(gene_data, chromosome, both_strands)
        if both_strands:
            _add_genes_to_G(G, genes, gene_data)
        else:
            for genes_per_strand in genes:
                _add_genes_to_G(G, genes_per_strand, gene_data)

    return G


def _get_genes_on_chromosome(gene_data, chromosome, both_strands):
    """Determines the genes located on the specified chromosome.

    Returns either two lists of genes (one for each strand of 'chromosome') if
    both_strands is set to False, or one list of genes (for the genes located on
    both strands of 'chromosome') if both_strands is set to True.

    :param gene_data: dict storing chromosome name, strand, and gene position,
        for every gene in the dict
    :param chromosome: name of a given chromosome
    :param both_strands: bool specifying whether genes on both strands of a
        given chromosome should be considered neighbors (if set to True) or only
        genes on a single strand should be considered neighbors (if set to
        False)
    :return: tuple of lists of genes (if both_strands is False) or list of genes
        (if both_strands is True)
    """
    genes = list()

    for gene in gene_data:
        if gene_data[gene]['chr'] == chromosome:
            genes.append(gene)

    genes = sorted(genes, key=lambda(g): gene_data[g]['pos'])

    if both_strands:
        return genes

    genes_fwd = [gene for gene in genes if gene_data[gene]['fwd']]
    genes_rev = [gene for gene in genes if gene not in genes_fwd]

    return genes_fwd, genes_rev


def _add_genes_to_G(G, genes, gene_data):
    """Adds the specified list fo genes to the graph G representing gene
    neighborhood.

    :param G: undirected graph representing gene neighborhood
    :param genes: list of genes to add to G
    :param gene_data: dict storing chromosome name, strand, and gene position,
        for every gene in the dict
    """
    for i in range(1, len(genes)):
        G.add_edge(genes[i], genes[i-1])

    if len(genes) > 1:
        if len(gene_data[genes[-1]]['pos']) > 1:
            start = gene_data[genes[-1]]['pos'][0][0]
            end = gene_data[genes[-1]]['pos'][-1][1]
            if start > end:
                G.add_edge(genes[0], genes[-1])


def add_edges_to_G(G_init, reactions, delta_G=0):
    """Builds and returns an undirected graph by adding supplementary edges to
    G_init, according to the value of the delta_G parameter.

    :param G_init: initial undirected graph to which supplementary edges are
        added if delta_G is greater than 0
    :param reactions: dict of dicts containing information on the reactions in
        a metabolic pathway, needed to select the vertices in G_init to be kept
    :param delta_G: integer designating how many genes CoMetGeNe can skip; if
        non-zero, edges are added to G_init between all pairs of vertices that
        are connected by a shortest path whose length is less than or equal to
        delta_G + 1
    :return: the undirected graph resulting from G_init with the newly added
        edges
    """
    G_init_nodes = G_init.nodes()
    involved = set()
    for r_id in reactions:
        for gene in reactions[r_id]['enzyme']:
            if gene in G_init_nodes:
                involved.add(gene)
            
    G = nx.Graph(nx.subgraph(G_init, list(involved)))

    if delta_G < 0:
        delta_G = 0

    if delta_G > 0:
        edges = G_init.edges()

        # Newly added edges are labeled with the list of skipped vertices in
        # G_init.
        for src in involved:
            for dst in involved:
                if src != dst and (src, dst) not in edges and \
                   nx.has_path(G_init, src, dst):
                    path = nx.shortest_path(G_init, src, dst)
                    skipped = path[(path.index(src)+1):path.index(dst)]
                    if 0 < len(skipped) <= delta_G:
                        G.add_edge(src, dst)
                        G[src][dst]['skipped'] = skipped

    return G


def add_edges_to_D(D_init, delta_D=0):
    """Adds arcs between vertices of D_init that are connected by a shortest
    path of length at most delta_D + 1, and returns the new graph.

    :param D_init: initial directed graph to which supplementary arcs are added
        if delta_D is greater than 0
    :param delta_D: integer designating how many reactions CoMetGeNe can skip;
        if non-zero, arcs are added to D_init between all pairs of vertices
        that are connected by a shortest path whose length is less than or equal
        to delta_D + 1
    :return: the directed graph resulting from D_init with the newly added arcs
    """
    D = nx.DiGraph(D_init)

    if delta_D < 0:
        delta_D = 0

    # Newly added arcs are labeled with the list of skipped vertices in D_init.
    if delta_D > 0:
        spl_D = dict(nx.shortest_path_length(D_init))
        for src in spl_D.keys():
            for dst in spl_D[src]:
                if src != dst and (src, dst) not in D.edges():
                    if spl_D[src][dst] <= delta_D + 1:
                        path = nx.shortest_path(D_init, src, dst)
                        D.add_edge(src, dst)
                        D[src][dst]['skipped'] = \
                            path[(path.index(src)+1):path.index(dst)]

    return D


def get_HNet_D(reactions, compounds):
    """Representing a metabolic pathway as a bipartite graph and returns the
    graph resulting from the projection of the bipartite graph on its set of
    reactions.

    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :param compounds: dict storing compound information (obtained by parsing a
        KGML file)
    :return: directed graph representing a metabolic pathway, in which vertices
        are reaction identifiers and arcs designate shared substrates between
        reactions
    """
    D = nx.DiGraph()
    D.add_nodes_from(reactions, bipartite=0)
    D.add_nodes_from(compounds, bipartite=1)

    for r_id in reactions:
        for substrate in reactions[r_id]['substrate']:
            D.add_edge(substrate, r_id)
            if reactions[r_id]['reversible']:
                D.add_edge(r_id, substrate)
        for product in reactions[r_id]['product']:
            D.add_edge(r_id, product)
            if reactions[r_id]['reversible']:
                D.add_edge(product, r_id)

    D_reduced = nx.algorithms.bipartite.projected_graph(D, reactions.keys())

    return D_reduced


def get_HNet_G(reactions, D_init, G_tmp):
    """Returns an undirected graph with reactions for vertices, representing
    gene neighborhood with respect to the reactions that the gene products
    catalyze.

    Two vertices ri and rj in this graph are connected by an edge if at least
    one of the genes coding for enzymes involved in reaction ri is neighboring a
    gene encoding an enzyme involved in rj. For more details, see Model in the
    methods section of Zaharia et al. (2018).

    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :param D_init: directed graph representing a metabolic pathway, in which
        vertices are reaction identifiers and arcs designate shared substrates
        between reactions; D_init does not contain additional arcs (for skipped
        reactions)
    :param G_tmp: directed graph representing gene neighborhood, in which
        vertices are gene names and edges translate neighboring genes; G_tmp
        contains additional edges (for skipped genes)
    :return: an undirected graph with reactions for vertices, representing gene
        neighborhood in terms of reactions
    """
    clean_up_reactions(reactions, D_init, G_tmp)

    G = nx.Graph()
    G.add_nodes_from(reactions.keys())

    enzyme_dict = dict()
    for r_id in reactions:
        for enzyme_name in reactions[r_id]['enzyme']:
            if enzyme_name not in enzyme_dict:
                enzyme_dict[enzyme_name] = list()
            if r_id not in enzyme_dict[enzyme_name]:
                enzyme_dict[enzyme_name].append(r_id)

    for (enzyme1, enzyme2) in G_tmp.edges():
        if enzyme1 in enzyme_dict.keys() and enzyme2 in enzyme_dict.keys():
            for r_id1 in enzyme_dict[enzyme1]:
                for r_id2 in enzyme_dict[enzyme2]:
                    if r_id1 != r_id2:
                        G.add_edge(r_id1, r_id2)
                        # Preserve information on skipped vertices.
                        if 'skipped' in G_tmp[enzyme1][enzyme2]:
                            G[r_id1][r_id2]['skipped'] = \
                                G_tmp[enzyme1][enzyme2]['skipped']

    return G


def clean_up_reactions(reactions, D_init, G_tmp):
    """Modifies reactions and D_init by removing superfluous enzyme names and
    internal KGML reaction IDs.

    Enzyme names are removed from reactions if no genes with the same name exist
    in G_tmp.

    Internal KGML reaction IDs represent keys in 'reactions' and vertices in
    D_init, respectively. They are removed if no genes associated to the
    reactions they designated exist in G_tmp.

    The 'reactions' dict and the D_init graph are modified.

    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :param D_init: directed graph representing a metabolic pathway, in which
        vertices are reaction identifiers and arcs designate shared substrates
        between reactions; D_init does not contain additional arcs (for skipped
        reactions)
    :param G_tmp: directed graph representing gene neighborhood, in which
        vertices are gene names and edges translate neighboring genes; G_tmp
        contains additional edges (for skipped genes)
    """
    candidate_ids_to_remove = set()
    for r_id in D_init.nodes():
        if r_id in reactions.keys():
            for enzyme_name in reactions[r_id]['enzyme']:
                if enzyme_name not in G_tmp.nodes():
                    candidate_ids_to_remove.add(r_id)

    for r_id in candidate_ids_to_remove:
        remove = True
        names_to_remove = set()
        for enzyme_name in reactions[r_id]['enzyme']:
            if enzyme_name in G_tmp.nodes():
                remove = False
                break
            else:
                names_to_remove.add(enzyme_name)
        if remove:  # No enzyme associated to r_id is found in G_tmp.
            del reactions[r_id]
            D_init.remove_node(r_id)
        else:  # Remove just enzyme names not found in G_tmp.
            for enzyme_name in names_to_remove:
                reactions[r_id]['enzyme'].remove(enzyme_name)
