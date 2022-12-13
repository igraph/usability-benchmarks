import igraph as ig

def num_incident_edges(g : ig.Graph) -> list[int]:
    n = len(g.vs)
    res = [0] * n
    for e in g.es:
        res[e.source] += 1
        if e.target != e.source:
            res[e.target] += 1
    return res

g = ig.Graph.Erdos_Renyi(n=6, m=12, directed=False, loops=True)
ig.plot(g, target='./solutions/python-igraph/number_of_incident_edges/graph.png', vertex_label=num_incident_edges(g), vertex_color='white')