import networkx as nx
import numpy as np
from scipy.optimize import linprog


def num_to_array(i, n):
    """Convert integer i to binary array s with length n"""
    """For example, s=1000...0 is i=1"""
    s = np.zeros(n, dtype=int)
    for j in range(0, n):
        s[j] = i % 2
        i = i / 2
    return s


def fS(G):
    """Generate f(S)"""
    n = G.number_of_nodes()
    f = np.zeros(2**n, dtype=int)
    for i in range(0, 2**n):
        s = num_to_array(i, n)
        for u in range(0, n):
            if (s[u] == 1) and (G.nodes[u]["worker"] == True):
                for v in range(0, n):
                    if (s[v] == 0) and (G.nodes[v]["worker"] == True):
                        f[i] = 1
    return f


def check_steiner_tree(F):
    """Check if F is a steiner tree"""
    n = F.number_of_nodes()
    for u in range(0, n):
        for v in range(u + 1, n):
            if (F.nodes[u]["worker"]) and (F.nodes[v]["worker"]):
                return nx.has_path(F, u, v)
    return True


def solve_lp(G, f):
    e = G.number_of_edges()
    n = G.number_of_nodes()

    c = np.zeros(e)
    for edge in G.edges:
        c[G.edges[edge]["index"]] = G.edges[edge]["weight"]

    A = np.zeros((2**n + n, e))
    edge_list = list(G.edges)

    for i in range(0, 2**n):
        s = num_to_array(i, n)
        for j in range(0, e):
            u, v = edge_list[j]
            if s[u] + s[v] == 1:
                A[i, j] = -1

    for i in range(2**n, 2**n + n):
        incident_edges = G.edges(i - 2**n)
        j = 0
        for edge in G.edges:
            if edge in incident_edges:
                A[i, j] = 1
            j = j + 1

    b = -f
    for node in G.nodes:
        b = np.append(b, G.nodes[node]["degree_constraint"])

    res = linprog(c, A, b)

    breakpoint()

    return res.x


def remove_degree_constraint(G):
    for node in G.nodes:
        if (G.degree(node) <= G.nodes[node]["degree_constraint"] + 4):
            G.nodes[node]["degree_constraint"] = None


def pick_one_edge(G, F, x):
    i = 0
    for edge in G.edges:
        if np.isclose(x[i], 1, rtol=1e-6):
            F.add_edge(edge[0], edge[1])
            G.remove_edge(edge[0], edge[1])
            u = edge[0]
            v = edge[1]
            if G.nodes[u]["degree_constraint"] != None:
                G.nodes[u]["degree_constraint"] -= 1
            if G.nodes[v]["degree_constraint"] != None:
                G.nodes[v]["degree_constraint"] -= 1
        i += 1


def rounding(G, F, x):
    i = 0
    for edge in G.edges:
        if x[i] >= 0.5:
            u = edge[0]
            v = edge[1]
            if (G.nodes[u]["degree_constraint"] == None) and (
                G.nodes[v]["degree_constraint"] == None):
                F.add_edge(edge[0], edge[1])
                G.remove_edge(edge[0], edge[1])
        i += 1


def f_update(f, F):
    n = F.number_of_nodes()
    for i in range(0, 2**n):
        if f[i] == 1:
            s = num_to_array(i, n)
            for u in range(0, n):
                if s[u] == 1:
                    for v in range(0, n):
                        if (s[v] == 0) and (v in F.neighbors(u)):
                            f[i] = 0


def dcstp(G):
    # add index to G's edges
    nx.set_edge_attributes(G, 0, "index")
    index = 0
    for edge in G.edges:
        G.edges[edge]["index"] = index
        index += 1

    F = nx.Graph()
    F.add_nodes_from(G.nodes.data())

    # find f(S)
    f = fS(G)

    while check_steiner_tree(F) == False:
        x = solve_lp(G, f)
        remove_degree_constraint(G)
        pick_one_edge(G, F, x)
        rounding(G, F, x)
        f_update(f, F)
    
    breakpoint()
    
    return F
