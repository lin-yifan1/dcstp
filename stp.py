import networkx as nx
import numpy as np
from dcstp import dcstp
import pdb

G = nx.Graph()
G.add_edges_from(np.genfromtxt("edge_list.txt", dtype=int), weight=2)

worker_list = np.atleast_1d(np.genfromtxt("worker_list.txt", dtype=int))

# add worker attribute to graph G
nx.set_node_attributes(G, False, "worker")
for worker in worker_list:
    G.nodes[worker]["worker"] = True

# node index starts from 0
# add degree_constraint attribute to graph G
n = G.number_of_nodes()
bv = int(input("Input degree constraint bv: "))
nx.set_node_attributes(G, bv, "degree_constraint")
for node in G.nodes:
    if G.nodes[node]['worker'] == True:
        G.nodes[node]['degree_constraint'] = 1


# setting edge weights to 1 or 2
for worker in worker_list:
    incident_edges = G.edges(worker)
    for edge in incident_edges:
        G.edges[edge]["weight"] = 1
    
F = dcstp(G)
