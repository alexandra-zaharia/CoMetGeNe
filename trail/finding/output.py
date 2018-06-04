"""
These utilities handle CoMetGeNe (trail finding) output.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

import networkx as nx
from trail import span, get_rn_trail


def output_trail(kgml, trail, net_inst, dev_out):
    """Outputs a given CoMetGeNe trail using internal KGML reaction IDs, R
    numbers, names of genes involved in the reactions, and EC numbers.

    Skipped reactions and/or genes are also displayed in the trail output.

    :param kgml: filename for metabolic pathway in KGML format
    :param trail: CoMetGeNe trail
    :param net_inst: NetworkBuilder object storing CoMetGeNe information and
        graphs
    :param dev_out: device for trail output (file or stdout)
    """
    dev_out.write(kgml + ':\n')

    found = kgml + \
        ": Found a trail of span %d " % \
        span(get_rn_trail(trail, net_inst.reactions))
    skipped_G = has_skipped_vertices_G(trail, net_inst.reactions, net_inst.G)
    skipped_D = has_skipped_vertices_D(trail, net_inst.D)

    if (net_inst.delta_G == 0 and net_inst.delta_D == 0) or \
            (not skipped_G and not skipped_D):
        found += 'without skipping any vertex'
    else:
        found += 'containing skipped vertices'

    dev_out.write(found + '\n')

    format_id = kgml + ': ' + format_trail(trail, net_inst.D)
    dev_out.write(format_id + '\n')
    format_rn = kgml + ': ' + format_trail(
        trail, net_inst.D, 'rn', net_inst.reactions)
    dev_out.write(format_rn + '\n')
    format_en = kgml + ': ' + format_trail(
        trail, net_inst.D, 'en', net_inst.reactions)
    dev_out.write(format_en + '\n')
    format_ec = kgml + ': ' + format_trail(
        trail, net_inst.D, 'ec', net_inst.reactions)
    dev_out.write(format_ec + '\n')

    if skipped_G:
        skipped = kgml + ': Skipped genes: '
        skipped += format_skipped_G(
            trail, net_inst.reactions, net_inst.G)
        dev_out.write(skipped + '\n')


def format_trail_helper_id(trail, D):
    """Formats a given CoMetGeNe trail using internal KGML reaction IDs.

    Skipped reactions are equally accounted for.

    :param trail: CoMetGeNe trail
    :param D: directed graph that served as input to the HNet algorithm,
        representing the reaction network in which trail was found
    :return: textual representation of the trail using internal KGML reaction
        IDs
    """
    string = ''

    for i in range(len(trail)):
        if len(string) > 0:
            string += ' -> '

        skipped = False
        if i > 0:
            if 'skipped' in D.edge[trail[i - 1]][trail[i]]:
                skipped = True

        if skipped:
            string += '['
            skipped_vertices = D.edge[trail[i - 1]][trail[i]]['skipped']
            for j in range(len(skipped_vertices)):
                string += skipped_vertices[j]
                if j < len(skipped_vertices) - 1:
                    string += ' -> '
            string += '] -> '

        string += str(trail[i])

    return string


def format_trail_helper_rn(trail, D, reactions):
    """Formats a given CoMetGeNe trail using R numbers.

    Skipped reactions are equally accounted for.

    :param trail: CoMetGeNe trail
    :param D: directed graph that served as input to the HNet algorithm,
        representing the reaction network in which trail was found
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :return: textual representation of the trail using R numbers
    """
    string = ''

    for i in range(len(trail)):
        if len(string) > 0:
            string += ' -> '

        skipped = False
        if i > 0:
            if 'skipped' in D.edge[trail[i - 1]][trail[i]]:
                skipped = True

        if skipped:
            string += '['
            skipped_vertices = D.edge[trail[i - 1]][trail[i]]['skipped']
            for j in range(len(skipped_vertices)):
                r_skipped_list = [
                    r.replace('rn:', '') for r in
                    reactions[skipped_vertices[j]]['reaction']
                ]
                r_skipped_path = ', '.join(r_skipped_list)
                if len(reactions[skipped_vertices[j]]['reaction']) > 1:
                    string += '{' + r_skipped_path + '}'
                else:
                    string += r_skipped_path
                if j < len(skipped_vertices) - 1:
                    string += ' -> '
            string += '] -> '

        r_list = [
            r.replace('rn:', '') for r in reactions[trail[i]]['reaction']
        ]
        r_path = ', '.join(r_list)
        if len(reactions[trail[i]]['reaction']) > 1:
            string += '{' + r_path + '}'
        else:
            string += r_path

    return string


def format_trail_helper_en(trail, D, reactions):
    """Formats a given CoMetGeNe trail using names of genes involved in every
    reaction in the trail.

    Skipped reactions are equally accounted for by showing which genes
    correspond to the skipped reactions.

    :param trail: CoMetGeNe trail
    :param D: directed graph that served as input to the HNet algorithm,
        representing the reaction network in which trail was found
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :return: textual representation of the trail using gene names of the genes
        involved in every reaction in the trail
    """
    string = ''

    for i in range(len(trail)):
        if len(string) > 0:
            string += ' -> '

        skipped = False
        if i > 0:
            if 'skipped' in D.edge[trail[i - 1]][trail[i]]:
                skipped = True

        if skipped:
            string += '['
            skipped_vertices = D.edge[trail[i - 1]][trail[i]]['skipped']
            for j in range(len(skipped_vertices)):
                e_skipped_path = \
                    ', '.join(reactions[skipped_vertices[j]]['enzyme'])
                if len(reactions[skipped_vertices[j]]['enzyme']) > 1:
                    string += '{' + e_skipped_path + '}'
                else:
                    string += e_skipped_path
                if j < len(skipped_vertices) - 1:
                    string += ' -> '
            string += '] -> '

        e_path = ', '.join(reactions[trail[i]]['enzyme'])
        if len(reactions[trail[i]]['enzyme']) > 1:
            string += '{' + e_path + '}'
        else:
            string += e_path

    return string


def format_trail_helper_ec(trail, D, reactions):
    """Formats a given CoMetGeNe trail using EC numbers associated to every
    reaction in the trail

    Skipped reactions are equally accounted for by showing which EC numbers
    correspond to the skipped reactions.

    :param trail: CoMetGeNe trail
    :param D: directed graph that served as input to the HNet algorithm,
        representing the reaction network in which trail was found
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :return: textual representation of the trail using EC numbers associated to
        every reaction in the trail
    """
    string = ''

    for i in range(len(trail)):
        if len(string) > 0:
            string += ' -> '

        skipped = False
        if i > 0:
            if 'skipped' in D.edge[trail[i - 1]][trail[i]]:
                skipped = True

        if skipped:
            string += '['
            skipped_vertices = D.edge[trail[i - 1]][trail[i]]['skipped']
            for j in range(len(skipped_vertices)):
                if len(reactions[skipped_vertices[j]]['ec']) > 0:
                    e_skipped_path = \
                        ', '.join(reactions[skipped_vertices[j]]['ec'])
                else:
                    e_skipped_path = '?'
                if len(reactions[skipped_vertices[j]]['ec']) > 1:
                    string += '{' + e_skipped_path + '}'
                else:
                    string += e_skipped_path
                if j < len(skipped_vertices) - 1:
                    string += ' -> '
            string += '] -> '

        if len(reactions[trail[i]]['ec']) > 0:
            e_path = ', '.join(reactions[trail[i]]['ec'])
        else:
            e_path = '?'
        if len(reactions[trail[i]]['ec']) > 1:
            string += '{' + e_path + '}'
        else:
            string += e_path

    return string


def format_trail(trail, D, formatting_type='id', reactions=None):
    """Returns a string representing the input CoMetGeNe trail, formatted
    according to the specified formatting_type.

    Possible values for formatting_type:
        'id' for internal KGML reaction IDs
        'rn' for R numbers
        'en' for gene names
        'ec' for EC numbers

    :param trail: CoMetGeNe trail
    :param D: directed graph that served as input to the HNet algorithm,
        representing the reaction network in which trail was found
    :param formatting_type: formatting to apply to the trail (possible values:
        'id', 'rn', 'en', 'ec')
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :return: textual representation of the trail in D, formatted as v1 -> ...
        -> vN, where vi is a numerical identifier (if formatting_type is 'id'),
        a (list of) KEGG reaction number(s) ('rn'), a (list of) gene name(s)
        ('en'), or a (list of) EC number(s)
    """
    if reactions is None:
        reactions = dict()
    if formatting_type not in ['id', 'rn', 'en', 'ec']:
        formatting_type = 'id'

    if formatting_type == 'id':
        return format_trail_helper_id(trail, D)
    elif formatting_type == 'rn':
        return format_trail_helper_rn(trail, D, reactions)
    elif formatting_type == 'en':
        return format_trail_helper_en(trail, D, reactions)
    else:
        return format_trail_helper_ec(trail, D, reactions)


def format_skipped_G(trail, reactions, G):
    """Returns a textual representation of the genes skipped to obtain the given
    CoMetGeNe trail.

    :param trail: CoMetGeNe trail
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :param G: undirected graph built on the same vertex set as the metabolic
        pathway (see Model in the methods section of Zaharia et al., 2018)
    :return: string of genes that were skipped in the gene neighborhood network
        in order to obtain the CoMetGeNe trail
    """
    skipped = set()

    G_sub = nx.subgraph(G, trail)
    nodes = G_sub.nodes()

    involved = set()  # genes involved in this pathway
    for r_id in nodes:
        for gene in reactions[r_id]['enzyme']:
            involved.add(gene)

    for v1 in G_sub.edge:
        for v2 in G_sub.edge[v1]:
            if 'skipped' in G_sub.edge[v1][v2]:
                for vertex in G_sub.edge[v1][v2]['skipped']:
                    if vertex not in involved:
                        skipped.add(vertex)

    return ', '.join(str(skipped_vertex) for skipped_vertex in skipped)


def has_skipped_vertices_G(trail, reactions, G):
    """Determines whether the given CoMetGene trail skips any genes.

    :param trail: CoMetGene trail
    :param reactions: dict of dicts storing reaction information (obtained by
        parsing a KGML file)
    :param G: undirected graph built on the same vertex set as the metabolic
        pathway (see Model in the methods section of Zaharia et al., 2018)
    :return: True if the CoMetGeNe trail was obtained skipping at least one
        gene, False otherwise
    """
    involved = set()  # genes involved in this trail
    for r_id in trail:
        for gene in reactions[r_id]['enzyme']:
            involved.add(gene)

    G_sub = nx.subgraph(G, trail)
    for v1 in G_sub.edge:
        for v2 in G_sub.edge[v1]:
            if 'skipped' in G_sub.edge[v1][v2]:
                for vertex in G_sub.edge[v1][v2]['skipped']:
                    if vertex not in involved:
                        return True

    return False


def has_skipped_vertices_D(trail, D):
    """Determines whether the given CoMetGene trail skips any reactions.

    :param trail: CoMetGene trail
    :param D: directed graph that served as input to the HNet algorithm,
        representing the reaction network in which trail was found
    :return: True if the CoMetGeNe trail was obtained skipping at least one
        reaction, False otherwise
    """
    for i in range(1, len(trail)):
        if 'skipped' in D.edge[trail[i-1]][trail[i]]:
            return True

    return False
