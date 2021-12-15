from igraph import Graph
import itertools


def k_cliques(g, k):
    """Returns a list of all cliques of size k in Graph g. Each clique is given as a set of vertex ids in g.
    """
    """
    edges = g.get_edgelist()
    cliques = [{i, j} for i, j in edges if i != j]
    for c in range(2, k):
        cliques_n = set()
        for i, j in itertools.combinations(cliques, 2):
            x = i ^ j
            if len(x) == 2 and (tuple(x) in edges or tuple(x)[::-1] in edges):
                cliques_n.add(tuple(i | j))
            
        cliques = list(map(set, cliques_n))
    
    return cliques
    """
    return g.cliques(k, k)


if __name__ == "__main__":
    g = Graph()
    g.add_vertices(7)
    g.add_edges([(0, 1), (1, 2), (1, 3), (1, 5), (2, 3), (2, 4), (2, 6), (3, 4), (3, 6), (4, 5), (4, 6)])
    print(k_cliques(g, 3))
