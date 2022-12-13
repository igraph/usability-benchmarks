import igraph as ig
import random

g = ig.Graph.Erdos_Renyi(n=20, p=0.4, directed=False, loops=False)
v = 5
res = g.distances(v, g.vs, weights=None, mode='all')
g.vs["dist_from_5"] = res[0]
ig.plot(g, target='./solutions/python-igraph/vertex_attributes/graph.png', vertex_label=g.vs["dist_from_5"], vertex_color='white')





