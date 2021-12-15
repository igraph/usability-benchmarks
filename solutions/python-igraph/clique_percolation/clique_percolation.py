from igraph import Graph

EXAMPLE_K = 3

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

def k_communities(g, k, min_common_vertices):
    """
    Given a graph and size k, return a list of communities found through clique
    percolation with at least min_common_vertices shared.
    A community is the connected component of cliques where two cliques overlap
    in at least k-1 vertices.
    Each community is a set of vertices in that community

    Approach: generate a new clique graph where each k-clique is a vertex, and 
    connect vertices if the cliques share min_common_vertices
    """
    if (min_common_vertices > k):
        raise ValueError("Error: min_common_vertices=%d must be greater than k=%d" % (min_common_vertices, k))

    # Each clique can be identified by it's index in 'cliques'
    cliques = [set(x) for x in k_cliques(g, k)]
    num_cliques = len(cliques)
    # Graph storing each clique as a vertex
    clique_graph = Graph()
    clique_graph.add_vertices(num_cliques)
    cg_edges = []
    # List of communities as a list of list of clique indices
    clique_communities = []
    # List of communities as a list of list of set of vertices
    communities = []

    # Get connected cliques
    for c1 in range(num_cliques):
        for c2 in range(c1, num_cliques):
            # If the two cliques share at least min_common_vertices,
            # they are part of the same community
            if c1 != c2 and len(cliques[c1].intersection(cliques[c2])) >= min_common_vertices:
                cg_edges.append((c1, c2))
    
    # Add edges for clique graph
    clique_graph.add_edges(cg_edges)

    # Merge cliques accordingly into communities
    to_merge = list(range(num_cliques))
    while (len(to_merge) > 0):
        # Get first unmerged vertex
        idx = to_merge[0]
        # Get all connected vertices
        vertices = clique_graph.subcomponent(idx, mode="all")
        # Remove all connected vertices
        for v in vertices:
            to_merge.remove(v)
        clique_communities.append(vertices)

    # Generate the list of clique communities as a list of set of vertices
    for cc in clique_communities:
        vertices = set()
        for c in cc:
            clique = cliques[c]
            vertices.update(clique)
        communities.append(vertices)

    return communities

if __name__ == "__main__":
    # Load in Zachary's karate club network
    print("For the Zachary Karate Club graph\n")
    g = Graph.Famous("Zachary")

    # Calculate cliques and communities
    print("Cliques of size k=%d" % EXAMPLE_K)
    print(k_cliques(g, EXAMPLE_K))
    print("")
    print("Communities from cliques of size k=%d, with %d shared vertices" % (EXAMPLE_K, EXAMPLE_K-1))
    print(k_communities(g, EXAMPLE_K, EXAMPLE_K-1))
