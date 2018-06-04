import networkx as nx


class AccessPoints(object):
    """Computes entry and exit points for every SCC in the line graph."""
    def __init__(self, LD):
        """Initializes the internal data structure needed to represent entry
        and exit points to and from any SCC in the line graph LD.

        :param LD: line graph for which entry and exit points will be determined
        """
        self.points = dict()  # access points to/from SCCs in a line graph

        C = nx.condensation(LD)
        for vertex in C.nodes():
            self.points[vertex] = dict()
            self.points[vertex]['entry'] = dict()
            self.points[vertex]['exit'] = dict()

        # For each arc (u, v) in C, determine the exit points of the SCC in the
        # line graph corresponding to vertex u, as well as the entry points to
        # the SCC corresponding to vertex v.
        for (u, v) in C.edges():
            U = C.node[u]['members']
            V = C.node[v]['members']

            # vertex_u and vertex_v are vertices in the line graph, so they
            # should be tuples (i.e. arcs of the initial directed graph).
            for vertex_u in U:
                assert type(vertex_u) is tuple
                for vertex_v in V:
                    assert type(vertex_v) is tuple

                    if vertex_u[1] == vertex_v[0]:
                        # vertex_u is an exit point in U towards V.
                        if v not in self.points[u]['exit']:
                            self.points[u]['exit'][v] = set()
                        self.points[u]['exit'][v].add(vertex_u)

                        # vertex_v is an entry point in V from U.
                        if u not in self.points[v]['entry']:
                            self.points[v]['entry'][u] = set()
                        self.points[v]['entry'][u].add(vertex_v)

    def entry_points(self, v_c, u_c):
        """Returns a list of entry points in the SCC corresponding to vertex v_c
        in the condensation graph when coming from the SCC corresponding to
        vertex u_c.

        :param v_c: vertex in the condensation graph (the "to" vertex)
        :param u_c: vertex in the condensation graph (the "from" vertex)
        :return: list of vertices in the line graph corresponding to entry
            points for the SCC whose condensation is v_c when coming from the
            SCC whose condensation is u_c
        """
        if u_c not in self.points[v_c]['entry']:
            return list()
        return self.points[v_c]['entry'][u_c]

    def exit_points(self, v_c, w_c):
        """Returns a list of exit points in the SCC corresponding to vertex v_c
        in the condensation graph when heading toward the SCC corresponding to
        vertex w_c.

        :param v_c: vertex in the condensation graph (the "from" vertex)
        :param w_c: vertex in the condensation graph (the "to" vertex)
        :return: list of vertices in the line graph corresponding to exit points
            for the SCC whose condensation is v_c when heading toward the SCC
            whose condensation is w_c
        """
        if w_c not in self.points[v_c]['exit']:
            return list()
        return self.points[v_c]['exit'][w_c]

    def entry_point_from(self, v_c, v_ld):
        """Returns the list of vertices in the condensation graph corresponding
        to SCCs in the line graph having outgoing edges that end in v_ld (i.e.
        v_ld is an entry point in the SCC corresponding to v_c with respect to
        the SCCs in the list returned by this method).

        :param v_c: vertex in the condensation graph
        :param v_ld: vertex in the line graph, more specifically in the SCC
            whose condensation is v_c
        :return: vertices in the condensation graph acting as starting points
            for outgoing edges that end in v_ld
        """
        from_sccs = list()
        for scc in self.points[v_c]['entry']:
            if v_ld in self.points[v_c]['entry'][scc]:
                from_sccs.append(scc)
        return from_sccs

    def exit_point_to(self, v_c, v_ld):
        """Returns the list of vertices in the condensation graph corresponding
        to SCCs in the line graph having incoming edges that start in v_ld (i.e.
        v_ld is an exit point in the SCC corresponding to v_c with respect to
        the SCCs in the list returned by this method).

        :param v_c: vertex in the condensation graph
        :param v_ld: vertex in the line graph, more specifically in the SCC
            whose condensation is v_c
        :return: vertices in the condensation graph acting as ending points for
            incoming edges that start in v_ld
        """
        to_sccs = list()
        for scc in self.points[v_c]['exit']:
            if v_ld in self.points[v_c]['exit'][scc]:
                to_sccs.append(scc)
        return to_sccs

    def entries_from_sccs(self, v_c):
        """Returns the list of vertices in the condensation graph representing
        predecessor SCCs to the SCC corresponding to v_c.

        :param v_c: vertex in the condensation graph
        :return: list of vertices in the condensation graph acting as
            predecessor SCCs for the SCC corresponding to v_c
        """
        return list(self.points[v_c]['entry'])

    def exits_to_sccs(self, v_c):
        """Returns the list of vertices in the condensation graph representing
        successor SCCs to the SCC corresponding to v_c.

        :param v_c: vertex in the condensation graph
        :return: list of vertices in the condensation graph acting as successor
            SCCs for the SCC corresponding to v_c
        """
        return list(self.points[v_c]['exit'])

    def has_entries(self, v_c):
        """Determines whether the SCC corresponding to vertex v_c in the
        condensation graph has at least one entry point.

        :param v_c: vertex in the condensation graph
        :return: True if the SCC corresponding to v_c has at least one entry
            point, False otherwise
        """
        return len(self.points[v_c]['entry']) > 0

    def has_exits(self, v_c):
        """Determines whether the SCC corresponding to vertex v_c in the
        condensation graph has at least one exit point.

        :param v_c: vertex in the condensation graph
        :return: True if the SCC corresponding to v_c has at least one exit
            point, False otherwise
        """
        return len(self.points[v_c]['exit']) > 0
