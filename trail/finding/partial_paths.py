"""
These utilities handle path finding in the line graph.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import networkx as nx

from trail import span, get_rn_trail, get_corresponding_trail


def search_cache(C, cache):
    """Determines whether the graph C exists in the cache.

    :param C: directed acyclic graph representing the condensation graph of a
        directed line graph
    :param cache: dict storing best partial paths for every strongly connected
        component (SCC) of the line graph whose condensation graph is C
    :return: the graph in the cache dict identical to C if it exists, None
        otherwise
    """
    for cond in cache.keys():
        if sorted(C.nodes()) == sorted(cond.nodes()):
            skip = False  # Skip graph 'cond'?
            for node in C.nodes():
                if C.node[node]['members'] != cond.node[node]['members']:
                    skip = True
                    break
            if not skip:  # Every SCC in C exists in cond (same members).
                return cond

    return None


def partial_paths(LD, access, cache, reactions):
    """Computes best partial paths for every SCC of LD, the line graph of the
    directed graph serving as input for the HNet algorithm. For more details see
    Path finding in the line graph in the methods section of Zaharia et al.
    (2018).

    The cache dict is modified, storing best partial paths in every SCC of LD
    between all possible pairs of entry and exit points from predecessor SCCs to
    successor SCCs.

    :param LD: line graph of the directed graph serving as input for the HNet
        algorithm
    :param access: AccessPoints object representing entry and exit points for
        every SCC of LD
    :param cache: dict storing best partial paths for every strongly connected
        component (SCC) of LD
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    """
    C = nx.condensation(LD)
    cond = search_cache(C, cache)
    if cond is None:
        cache[C] = dict()
        cond = C
    for x in C.nodes():
        if x not in cache[cond]:
            cache[cond][x] = dict()

    if len(C.nodes()) == 1:
        c_node = C.nodes()[0]
        for src in LD.nodes():
            for dst in LD.nodes():
                if src == dst:
                    evaluate_path(
                        [src], None, None, cache[cond][c_node], reactions)
                else:
                    for path in nx.all_simple_paths(LD, src, dst):
                        evaluate_path(
                            path, None, None, cache[cond][c_node], reactions)
    else:
        for x in C.nodes():
            scc_x = nx.subgraph(LD, C.node[x]['members'])
            partial_paths_helper(x, scc_x, access, cache[cond][x], reactions)


def partial_paths_helper(x, scc_x, access, cache_scc, reactions):
    """Determines best partial paths in the SCC scc_x of the line graph.

    The cache dict is modified, storing best partial paths in every SCC of LD
    between all possible pairs of entry and exit points from predecessor SCCs to
    successor SCCs.

    :param x: vertex in the condensation graph of a line graph
    :param scc_x: SCC in the line graph of a directed graph whose condensation
        is vertex x
    :param access: AccessPoints object representing entry and exit points for
        every SCC of the line graph
    :param cache_scc: cache entry corresponding to vertex x in the condensation
        graph of the line graph
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    """
    sccs_from = access.entries_from_sccs(x)
    sccs_to = access.exits_to_sccs(x)

    if len(sccs_from) == len(sccs_to) == 0:
        for src in scc_x.nodes():
            for dst in scc_x.nodes():
                if src == dst:
                    evaluate_path([src], None, None, cache_scc, reactions)
                else:
                    for path in nx.all_simple_paths(scc_x, src, dst):
                        evaluate_path(path, None, None, cache_scc, reactions)
    else:
        partial_paths_between_access_points(
            x, scc_x, access, cache_scc, reactions)


def partial_paths_between_access_points(x, scc_x, access, cache_scc, reactions):
    """Determines best partial paths in the SCC scc_x of the line graph, if
    scc_x has at least one predecessor SCC or at least one successor SCC.

    The cache dict is modified, storing best partial paths in every SCC of LD
    between all possible pairs of entry and exit points from predecessor SCCs to
    successor SCCs.

    :param x: vertex in the condensation graph of a line graph
    :param scc_x: SCC in the line graph of a directed graph whose condensation
        is vertex x
    :param access: AccessPoints object representing entry and exit points for
        every SCC of the line graph
    :param cache_scc: cache entry corresponding to vertex x in the condensation
        graph of the line graph
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    """
    sccs_from = access.entries_from_sccs(x)
    sccs_to = access.exits_to_sccs(x)

    # From all possible entry points to all possible exit points in scc_x
    for w in sccs_from:
        for y in sccs_to:
            for src in access.entry_points(x, w):
                for dst in access.exit_points(x, y):
                    if src == dst:
                        evaluate_path([src], w, y, cache_scc, reactions)
                    else:
                        for path in nx.all_simple_paths(scc_x, src, dst):
                            evaluate_path(path, w, y, cache_scc, reactions)

    # From all possible entry points to every vertex in scc_x
    for w in sccs_from:
        for src in access.entry_points(x, w):
            for dst in scc_x.nodes():
                if src == dst:
                    evaluate_path([src], w, None, cache_scc, reactions)
                else:
                    for path in nx.all_simple_paths(scc_x, src, dst):
                        evaluate_path(path, w, None, cache_scc, reactions)

    # From every vertex in scc_x to all possible exit points
    for src in scc_x.nodes():
        for y in sccs_to:
            for dst in access.exit_points(x, y):
                if src == dst:
                    evaluate_path([src], None, y, cache_scc, reactions)
                else:
                    for path in nx.all_simple_paths(scc_x, src, dst):
                        evaluate_path(path, None, y, cache_scc, reactions)


def evaluate_path(path, scc_from, scc_to, cache_scc, reactions):
    """Evaluates a given path in a SCC of the line graph in terms of span and
    length in order to determine whether it is the best partial path so far in
    the SCC, from the predecessor SCC corresponding to scc_from and toward the
    successor SCC corresponding to scc_to.

    The cache dict is modified, storing best partial paths in every SCC of LD
    between all possible pairs of entry and exit points from predecessor SCCs to
    successor SCCs.

    :param path: path in a SCC of the line graph
    :param scc_from: vertex in the condensation graph of the line graph
        corresponding to a SCC in the line graph with an outgoing arc ending in
        the first vertex of path (can be None in case there is no predecessor
        SCC)
    :param scc_to: vertex in the condensation graph of the line graph
        corresponding to a SCC in the line graph with an incoming arc starting
        in the last vertex of path (can be None in case there is no successor
        SCC)
    :param cache_scc: cache entry for the vertex in the condensation graph of
        the line graph corresponding to the SCC in which path is enumerated
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    """
    evaluate_path_helper(path, (None, None), (None, None), cache_scc, reactions)
    if scc_from is not None or scc_to is not None:
        pts_pair = (path[0], path[-1])
        evaluate_path_helper(
            path, (scc_from, scc_to), pts_pair, cache_scc, reactions)


def evaluate_path_helper(path, scc_pair, pts_pair, cache_scc, reactions):
    """Evaluates a given path in a SCC of the line graph in terms of span and
    length in order to determine whether it is the best partial path so far in
    the SCC, between SCCs designated by scc_pair and between vertices designated
    by pts_pair.

    The cache dict is modified, storing best partial paths in every SCC of LD
    between all possible pairs of entry and exit points from predecessor SCCs to
    successor SCCs.

    :param path: path in a SCC of the line graph
    :param scc_pair: tuple of vertices in the condensation graph of the line
        graph, designating the predecessor and successor SCCs of the SCC in
        which path is enumerated
    :param pts_pair: tuple of vertices in the line graph representing the first
        and last vertices of path
    :param cache_scc: cache entry for the vertex in the condensation graph of
        the line graph corresponding to the SCC in which path is enumerated
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    """
    if scc_pair not in cache_scc:
        cache_scc[scc_pair] = dict()
    if pts_pair not in cache_scc[scc_pair]:
        cache_scc[scc_pair][pts_pair] = path
    else:
        max_path = cache_scc[scc_pair][pts_pair]
        cache_scc[scc_pair][pts_pair] = best_path([max_path, path], reactions)


def find_paths(LD, cond_path, cache_ld):
    """Finds all paths in the line graph LD corresponding to a path cond_path
    in the condensation graph of LD.

    :param LD: line graph of the directed graph serving as input for the HNet
        algorithm
    :param cond_path: path in the condensation graph of LD
    :param cache_ld: cache entry for the vertex in the condensation graph of
        LD corresponding to the whole line graph LD
    :return: all paths in LD corresponding to the path cond_path in the
        condensation graph of LD
    """
    all_paths = list()
    find_paths_helper(LD, cond_path, 0, list(), all_paths, cache_ld)
    return all_paths


def find_paths_helper(LD, cond_path, index, partial, all_paths, cache_ld):
    """Finds all paths in the line graph LD corresponding to a path cond_path
    in the condensation graph of LD.

    The 'all_paths' list is modified, storing in the end all paths in LD that
    correspond to path cond_path in the condensation graph of LD

    :param LD: line graph of the directed graph serving as input for the HNet
        algorithm
    :param cond_path: path in the condensation graph of LD
    :param index: index of the current vertex examined in cond_path
    :param partial: partial path in LD corresponding to the first 'index'
        vertices of cond_path
    :param all_paths: list of all paths in LD that correspond to path cond_path
        in the condensation graph of LD
    :param cache_ld: cache entry for the vertex in the condensation graph of
        LD corresponding to the whole line graph LD
    """
    if index == len(cond_path):
        all_paths.append(partial)
    else:
        predecessor = cond_path[index - 1] if index > 0 else None
        successor = cond_path[index + 1] if index < len(cond_path) - 1 else None
        scc_pair = (predecessor, successor)

        if scc_pair in cache_ld[cond_path[index]]:
            for pts_pair in cache_ld[cond_path[index]][scc_pair]:
                if len(partial) == 0 or LD.has_edge(partial[-1], pts_pair[0]):
                    partial_concat = list(partial)
                    partial_concat.extend(
                        cache_ld[cond_path[index]][scc_pair][pts_pair])
                    find_paths_helper(
                        LD, cond_path, index + 1, partial_concat, all_paths,
                        cache_ld)


def best_path(candidates, reactions):
    """Determines and returns the best path among a list of candidate paths in
    the line graph, i.e. a path whose corresponding trail in the directed graph
    has maximum span, or maximum span and minimum length.

    :param candidates: list of candidate paths in the line graph
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :return: best path among the candidates (i.e. having a corresponding trail
        in the directed graph with maximum span, or with maximum span and
        minimum length)
    """
    if len(candidates) == 0:
        return list()

    best = list()

    for path in candidates:
        rn_trail = get_rn_trail(get_corresponding_trail(path), reactions)
        rn_best = get_rn_trail(get_corresponding_trail(best), reactions)
        if span(rn_trail) > span(rn_best):
            best = path
        elif span(rn_trail) == span(rn_best) and len(rn_trail) < len(rn_best):
            best = path

    return best
