"""
Utilities for cover set computation.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

from os.path import basename

import networkx as nx

from CoMetGeNeError import CoMetGeNeError


def path_in_graph(graph, path):
    """Determines whether the path is valid in the given undirected graph.

    :param graph: undirected graph
    :param path: list of vertices in the graph
    :return: True if the path is valid within the given graph, False otherwise
    """
    assert len(path) >= 2

    nodes = graph.nodes()
    for i in range(1, len(path)):
        if not (path[i] in nodes and
                path[i-1] in nodes and
                graph.has_edge(path[i-1], path[i])):
            return False

    return True


def path_in_cc(path, cc):
    """Determines whether all vertices of a given path belong to a given
    connected component.

    :param path:
    :param cc: list of vertices representing a connected component in a graph
    :return: True if the path vertices are found in the connected component,
        False otherwise
    """
    for node in path:
        if node not in cc:
            return False
    return True


def ccc(D, G, path):
    """Computes and returns the maximal common connected component (CCC)
    including vertices in 'path' of the two undirected graphs D and G.

    :param D: undirected graph with the same vertex set as G
    :param G: undirected graph with the same vertex set as D
    :param path: path in D
    :return: maximum CCC of D and G, including vertices in 'path'
    """
    assert type(D) is nx.Graph and type(G) is nx.Graph
    assert len(path) >= 2

    if not path_in_graph(D, path):
        return list()

    partitions = list()  # List of partitions of D and G, initially empty.
    ccc_helper(D, G, path, partitions)

    max_ccc = list()  # The CCC of D and G of maximum size, initially empty.
    for partition in partitions:
        if len(max_ccc) < len(partition):
            max_ccc = partition

    return max_ccc


def ccc_helper(D, G, path, partitions):
    """Computes all the common connected components (CCCs) of graphs D and G
    that include every vertex in 'path'. These CCCs are stored in the list
    'partitions'.

    The 'partitions' list is modified.

    :param D: directed graph with the same vertex set as G
    :param G: undirected graph with the same vertex set as D
    :param path: path in D
    :param partitions: initially empty list that will store CCCs of D and G that
        include every vertex in 'path'
    """
    ccd = list(nx.connected_components(D))
    ccg = list(nx.connected_components(G))

    for d in ccd:
        for g in ccg:
            if path_in_cc(path, d) and path_in_cc(path, g):
                inter = list(set(d) & set(g))
                D_ = nx.subgraph(D, inter)
                G_ = nx.subgraph(G, inter)
                if nx.is_connected(D_) and nx.is_connected(G_):
                    partitions.append(inter)
                else:
                    ccc_helper(D_, G_, path, partitions)


def get_bridges(D, G, path):
    """Computes and returns a list of bridges of 'path' in graph D with respect
    to G, where D and G are built on the same vertex set V.

    A bridge of 'path' in D with respect to G is a vertex r in V such that there
    is no common connected component of D*[V\{r}] and G[V\{r}] containing all
    the vertices in 'path'. D* stands for the undirected graph underlying D
    (obtained from D by removing arc orientation). H[X] is the subgraph of H
    induced by vertices in X.

    :param D: directed graph with the same vertex set as G
    :param G: undirected graph with the same vertex set as G
    :param path: path in D
    :return: list of bridges of 'path' in D w.r.t. G
    """
    assert sorted(D.nodes()) == sorted(G.nodes())
    bridges = list()

    for r in G.nodes():
        D_tmp = nx.Graph(D)
        G_tmp = nx.Graph(G)
        for graph in D_tmp, G_tmp:
            graph.remove_node(r)
        if len(ccc(D_tmp, G_tmp, path)) == 0:
            bridges.append(r)

    return bridges


def get_cover_set(D, G, path):
    """Computes and returns the cover set of 'path' in the directed graph D and
    in the undirected graph G as described in Fertin et al. (2012) Algorithms
    for subnetwork mining in heterogeneous networks.

    :param D: directed graph with the same vertex set as G
    :param G: undirected graph with the same vertex set as D
    :param path: path in D
    :return: the cover set of 'path' in D w.r.t. G
    """
    for vertex in path:  # ensure all vertices in 'path' exist in D
        if not D.has_node(vertex):
            raise CoMetGeNeError('hnet_node', basename(__file__), vertex)

    for i in range(1, len(path)):  # ensure 'path' is a valid path in D
        if not D.has_edge(path[i-1], path[i]):
            raise CoMetGeNeError(
                'hnet_edge', basename(__file__), path[i-1], path[i], path)

    init_set = set()
    for vtx_set in nx.ancestors(D, path[0]), path, nx.descendants(D, path[-1]):
        for vertex in vtx_set:
            init_set.add(vertex)

    D_S = nx.Graph(nx.subgraph(D, init_set))  # undirected graph underlying D[S]
    G_S = nx.subgraph(G, init_set)
    cover_set = ccc(D_S, G_S, path)
    stop = False
    while not stop and len(cover_set) > 0:
        cs_tmp = cover_set
        bridges = get_bridges(D, G, path)
        for r in bridges:
            bridge_cover = [r]
            for vertices in nx.ancestors(D, r), nx.descendants(D, r):
                for vertex in vertices:
                    bridge_cover.append(vertex)
            cs_tmp = list(set(cs_tmp) & set(bridge_cover))
        if sorted(cover_set) == sorted(cs_tmp):
            stop = True
        else:
            cover_set = cs_tmp
            D_S = nx.Graph(nx.subgraph(D, cover_set))
            G_S = nx.subgraph(G, cover_set)
            cover_set = ccc(D_S, G_S, path)

    return cover_set
