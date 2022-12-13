import igraph as ig

def mean_neighborhood_degree(g: ig.Graph, k: int):
    n = len(g.vs)
    res = [0] * n
    for v in g.vs:
        kneighborhood = []
        distances = g.distances(v, g.vs, weights=None, mode='all')[0]
        for i in range(n):
            if distances[i] == k:
                kneighborhood.append(i)
        for i in kneighborhood:
            res[v.index] += g.degree(i)
        res[v.index] /= len(kneighborhood)
        res[v.index] = round(res[v.index], 2)
    return res

g = ig.Graph.Famous("Zachary")
ig.plot(g, target='./solutions/python-igraph/mean_degree_neighborhood/graph1.png', vertex_size=40, vertex_label=mean_neighborhood_degree(g, 1), vertex_color=["white"]*len(g.vs), layout=g.layout_circle())
ig.plot(g, target='./solutions/python-igraph/mean_degree_neighborhood/graph2.png', vertex_size=40, vertex_label=mean_neighborhood_degree(g, 2), vertex_color=["white"]*len(g.vs), layout=g.layout_circle())
