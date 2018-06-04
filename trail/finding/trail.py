"""
Utilities for handling CoMetGeNe trails.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""


def span(lst):
    """Returns the number of unique vertices in the input list.

    :param lst: list of vertices in a graph or of tuples of vertices (i.e. arcs)
        in its line graph
    :return: number of unique vertices in the list
    """
    if len(lst) == 0:
        return 0
    if type(lst[0]) is not tuple:
        return len(set(lst))

    # lst is a list of tuples.
    vertex_set = set(list(lst[0]))
    for arc in lst[1:]:
        vertex_set.add(arc[1])
    return len(vertex_set)


def _get_rn_trail(trail, reactions):
    """Returns a textual representation of the CoMetGeNe trail using R numbers
    instead of the default internal KGML reaction IDs.

    :param trail: CoMetGeNe trail
    :param reactions:  dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :return: textual representation of R numbers representing the CoMetGeNe
        trail in terms of R numbers instead of internal KGML reaction IDs
    """
    rn_trail = list()

    for vertex in trail:
        if len(reactions[vertex]['reaction']) == 1:
            rn_trail.append(reactions[vertex]['reaction'][0])
        else:
            rn_trail.append(
                '{' + ', '.join(reactions[vertex]['reaction']) + '}')

    return rn_trail


def get_rn_trail(trail, reactions):
    """Given a CoMetGeNe trail in a directed graph (or a path in the line graph
    of the directed graph), returns the corresponding CoMetGeNe trail using R
    numbers instead of the default internal KGML reaction IDs.

    :param trail: CoMetGeNe trail or path in the line graph of the directed
        graph representing the reaction network for trail finding
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :return: textual representation of R numbers representing the CoMetGeNe
        trail (or the path in the line graph of the directed graph) in terms of
        R numbers instead of internal KGML reaction IDs
    """
    if len(trail) == 0:
        return list()

    if type(trail[0]) is not tuple:
        return _get_rn_trail(trail, reactions)
    else:  # 'trail' is a path in the line graph
        return _get_rn_trail(get_corresponding_trail(trail), reactions)


def get_corresponding_trail(path):
    """Given a path in the line graph of a directed graph, returns the trail in
    the directed graph corresponding to the path.

    :param path: path in a line graph of a directed graph
    :return: the trail in the directed graph corresponding to the path in the
        line graph
    """
    trail = list()

    if len(path) > 0:
        trail = list(path[0])
        for arc in path[1:]:
            trail.append(arc[1])

    return trail
