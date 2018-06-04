"""
Utilities to check whether CoMetGeNe trails returned by the HNet algorithm are
consistent within the original graphs (gene network and reaction network).

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""
import networkx as nx


def is_consistent(queue, G, reactions, trail):
    """Checks whether the trail also makes sense when between a metabolic
    pathway and the original gene neighborhood G.

    :param queue: multiprocessing queue for storing the result of the procedure
    :param G: undirected graph representing a genome, with genes
            for vertices (i.e. not reaction identifiers as are vertices in D).
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :param trail: CoMetGeNe trail
    :return: True if trail is consistent with respect to the gene neighborhood
        graph G, False otherwise
    """
    G_nodes = G.nodes()
    gene_trail = list()
    for vertex in trail:
        gene_vertex = list()
        for gene in reactions[vertex]['enzyme']:
            if gene in G_nodes:
                gene_vertex.append(gene)
        gene_trail.append(gene_vertex)

    decompositions = decompose(gene_trail)
    for decomp in decompositions:
        for X in nx.connected_components(nx.subgraph(G, decomp)):
            skip = False
            for vertex in trail:
                if not skip:
                    enzyme_names = reactions[vertex]['enzyme']
                    if len(set(X) & set(enzyme_names)) == 0:
                        skip = True
            if not skip:
                queue.put(True)
                return

    queue.put(False)


def decompose(gene_trail):
    """Performs all possible decompositions of a list of lists, taking exactly
    one item from each list.

    For example, if gene_trail is [['a', 'b'], ['c'], ['d','e']], all possible
    decompositions having exactly one item from each list are:
        ['a', 'c', 'd']
        ['a', 'c', 'e']
        ['b', 'c', 'd']
        ['b', 'c', 'e']

    :param gene_trail: list of lists where list items are identifiers of genes
        involved in reactions of a CoMetGeNe trail
    :return: all possible decompositions of gene_trail taking exactly one item
        from each list, in the form of a generator
    """
    return decompose_helper(gene_trail, dict())


def decompose_helper(gene_trail, selection):
    """Given a list of lists gene_trail, determines for every list in gene_trail
    the index of the item to be selected in that list such that a unique
    decomposition of gene_trail (consisting in exactly one item from each list)
    is obtained.

    For example, if gene_trail is [['a', 'b'], ['c'], ['d','e']], then:
        indexes [0, 0, 1] correspond to the decomposition ['a', 'c', 'e']
        indexes [1, 0, 1] correspond to the decomposition ['b', 'c', 'e']

    The 'selection' dict is modified.

    :param gene_trail: list of lists where list items are identifiers of genes
        involved in reactions of a CoMetGeNe trail
    :param selection: dict storing the currently selected item for every list in
        gene_trail
    :return: all possible decompositions of gene_trail taking exactly one item
        from each list, in the form of a generator
    """
    unique_genes = True

    for i in range(len(gene_trail)):
        if type(gene_trail[i]) is str:
            selection[i] = 0
        else:
            if i not in selection.keys():
                unique_genes = False
                for j in range(len(gene_trail[i])):
                    selection[i] = j
                    if i < len(gene_trail)-1:
                        for k in range(i+1, len(gene_trail)):
                            selection.pop(k, None)
                    for decomp in decompose_helper(gene_trail, selection):
                        yield decomp

    if unique_genes:
        yield do_decomposition(gene_trail, selection)


def do_decomposition(gene_trail, selection):
    """Given a list of lists gene_trail and indexes for every item to be
    selected in every list in gene_trail, returns a list representing the
    corresponding decomposition.

    For example, if gene_trail is [['a', 'b'], ['c'], ['d','e']] and the index
    for the first list (['a', 'b']) is 0, the index for the second list (['c'])
    is 0, and the index for the third list (['d', 'e']) is 1, then the
    corresponding decomposition of gene_trail is ['a', 'c', 'e'].

    :param gene_trail: list of lists where list items are identifiers of genes
        involved in reactions of a CoMetGeNe trail
    :param selection: dict storing the currently selected item for every list in
        gene_trail
    :return: list representing the decomposition of gene_trail given by indexes
        stored in 'selection'
    """
    for i in range(len(gene_trail)):
        assert i in selection.keys()
    for list_index in selection:
        assert 0 <= selection[list_index] < len(gene_trail[list_index])

    decomposition = list()

    for i in range(len(gene_trail)):
        decomposition.append(gene_trail[i][selection[i]])

    return decomposition
