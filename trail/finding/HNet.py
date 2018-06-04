"""
Implements the HNet algorithm as described in Zaharia et al. (2018).

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import networkx as nx

from AccessPoints import AccessPoints
from partial_paths import partial_paths, search_cache, best_path, find_paths
from trail import span, get_rn_trail, get_corresponding_trail
from cover_set import get_cover_set


def HNet_on_every_arc(queue, G, D, reactions):
    """Performs trail finding for graphs D and G built on the same vertex set,
    as described in Zaharia et al. (2018), for every arc in D.

    Trails of maximum span passing through arcs in D such that their vertex sets
    induce connected subgraphs in G are put in the multiprocessing queue
    received as argument.

    The queue argument is modified.

    :param queue: multiprocessing queue for storing all CoMetGeNe trails for
        graphs D and G
    :param G: undirected graph with the same vertex set as D
    :param D: directed graph with the same vertex set as G
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    """
    cache = dict()  # will store best partial paths for SCCs of L(D)
    all_trails = list()  # HNet results in descending order of trail span

    for arc in D.edges():
        trail = HNet(D, G, arc, cache, reactions)
        if len(trail) > 0 and trail not in all_trails:
            if can_add_trail(trail, all_trails, reactions):
                if len(all_trails) == 0:
                    all_trails.append(trail)
                else:
                    index = 0
                    while index < len(all_trails) and \
                            span(get_rn_trail(trail, reactions)) < \
                            span(get_rn_trail(all_trails[index], reactions)):
                        index += 1
                    all_trails.insert(index, trail)

    queue.put(all_trails)


def can_add_trail(candidate, results, reactions):
    """Determines whether a candidate trail can be added to CoMetGeNe results
    for the same input graphs.

    :param candidate: CoMetGeNe trail
    :param results: list of CoMetGeNe trails
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :return: True if the candidate trail has at least span 2 and if its vertices
        are not a permutation of vertices in another trail in 'results'
    """
    if span(get_rn_trail(candidate, reactions)) == 1:
        return False

    for trail in results:
        # Exclude candidate trail if the set of its vertices is a strict
        # permutation as that of another trail.
        if set(candidate) == set(trail):
            return False

    return True


def HNet(D, G, arc, cache, reactions):
    """Performs trail finding as described in Zaharia et al. (2018).

    :param D: directed graph with the same vertex set as G
    :param G: undirected graph with the same vertex set as D
    :param arc: arc in D
    :param cache: dict storing best partial paths for every SCC of the line
        graph of D
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :return: a trail of maximum span in D such that its vertex set induces a
        connected subgraph in G
    """
    assert type(D) is nx.DiGraph and type(G) is nx.Graph
    assert len(arc) == 2
    S = get_cover_set(D, G, list(arc))
    return max_span_trail(
        nx.subgraph(D, S), nx.subgraph(G, S), arc, cache, reactions)


def max_span_trail(D, G, input_arc, cache, reactions):
    """Given a directed graph D and an undirected graph G built on the same
    vertex set and in input arc in D, returns a trail of maximum span in D that
    passes through the input arc, such that its vertex set induces a connected
    subgraph in G.

    The cache dict is modified, storing best partial paths in every SCC of LD
    between all possible pairs of entry and exit points from predecessor SCCs to
    successor SCCs.

    :param D: input directed graph for the HNet algorithm, with the same vertex
        set as G and representing a metabolic pathway
    :param G: input undirected graph for the HNet algorithm, with the same
        vertex set as D and representing gene neighborhood (in terms of
        reactions; see Model in the methods section of Zaharia et al., 2018)
    :param input_arc: input arc in D on which the HNet algorithm is executed
    :param cache: dict storing best partial paths for every strongly connected
        component (SCC) of the line graph of D
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :return: a trail of maximum span in D such that its vertex set induces a
        connected subgraph in G
    """
    LD = nx.line_graph(D)
    if nx.number_of_nodes(LD) == 0:
        return list()

    access = AccessPoints(LD)  # access points for every SCC in the line graph
    partial_paths(LD, access, cache, reactions)

    return get_corresponding_trail(
        max_span_path(LD, G, input_arc, cache, reactions))


def max_span_path(LD, G, input_arc, cache, reactions):
    """Returns a path in the line graph LD of a directed graph whose
    corresponding trail in the directed graph has maximum span, passes through
    the input arc, and induces a connected subgraph in G.

    :param LD: line graph of the directed graph serving as input for the HNet
        algorithm
    :param G: input undirected graph for the HNet algorithm, with the same
        vertex set as the directed graph whose line graph is LD and representing
        gene neighborhood (in terms of reactions; see Model in the methods
        section of Zaharia et al., 2018)
    :param input_arc: input arc in the directed graph built on the same vertex
        set as G and serving as input for the HNet algorithm
    :param cache: dict storing best partial paths for every strongly connected
        component (SCC) of LD
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :return: a path in LD with a corresponding trail of maximum span in the
        directed graph that passes through input_arc and induces a connected
        subgraph in G
    """
    C = nx.condensation(LD)
    cond = search_cache(C, cache)
    assert cond is not None

    nodes = cond.nodes()
    vertex = None
    for cond_vertex in nodes:
        if input_arc in cond.node[cond_vertex]['members']:
            vertex = cond_vertex
            break
    assert vertex is not None

    path = list()
    for src in nodes:
        for dst in nodes:
            if src == dst == vertex:
                P = cache[cond][vertex][(None, None)][(None, None)]
                if input_arc in P and connected(G, P):
                    path = best_path([path, P], reactions)
            else:
                for Q in nx.all_simple_paths(cond, src, dst):
                    if vertex in Q:
                        for P in find_paths(LD, Q, cache[cond]):
                            if input_arc in P and connected(G, P):
                                path = best_path([path, P], reactions)

    return path


def connected(G, line_path):
    """Determines whether the vertex set of the trail corresponding to a path in
    the line graph induces a connected subgraph in G.

    :param G: input undirected graph for the HNet algorithm, with the same
        vertex set as the directed graph and representing gene neighborhood (in
        terms of reactions; see Model in the methods section of Zaharia et al.,
        2018)
    :param line_path: path in the line graph of the directed graph that is built
        on the same vertex set as G and that serves as input for the HNet
        algorithm
    :return: True if the subgraph induced in G by the vertex set of line_path
        is connected, False otherwise
    """
    subgraph = nx.subgraph(G, get_corresponding_trail(line_path))
    if subgraph.number_of_nodes() == 0:
        return False
    return nx.is_connected(subgraph)
