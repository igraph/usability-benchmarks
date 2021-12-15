import igraph as ig
import itertools
import matplotlib.pyplot as plt

EXAMPLE_K = 3

def k_cliques(g, k):
    """Returns a list of all cliques of size k in Graph g. Each clique is given as a set of vertex ids in g.
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
    cliques = list(map(set, k_cliques(g, k)))
    num_cliques = len(cliques)
    # Graph storing each clique as a vertex
    clique_graph = ig.Graph()
    clique_graph.add_vertices(num_cliques)
    # List of communities as a list of list of set of vertices
    communities = []

    # Get connected cliques
    # If the two cliques share at least min_common_vertices,
    # they are part of the same community
    cg_edges = [(c1, c2) for c1, c2 in itertools.combinations(range(num_cliques), 2) if len(cliques[c1].intersection(cliques[c2])) >= min_common_vertices]

    # Add edges for clique graph
    clique_graph.add_edges(cg_edges)

    # Merge cliques accordingly into communities
    cg_clusters = clique_graph.clusters()

    # Generate the list of clique communities as a list of set of vertices
    for g in cg_clusters:
        vertices = set()
        for c in g:
            vertices.update(cliques[c])
        communities.append(vertices)

    return communities

if __name__ == "__main__":
    # Load in Zachary's karate club network
    print("For the Zachary Karate Club graph\n")
    g = ig.Graph.Famous("Zachary")

    # Calculate cliques and communities
    communities = k_communities(g, EXAMPLE_K, EXAMPLE_K-1)
    print("Cliques of size k=%d" % EXAMPLE_K)
    print(k_cliques(g, EXAMPLE_K))
    print("")
    print("Communities from cliques of size k=%d, with %d shared vertices" % (EXAMPLE_K, EXAMPLE_K-1))
    print(communities)
    
    # plot the graph with a vertex cover for each community
    fig, ax = plt.subplots()
    ig.plot(ig.VertexCover(g, communities), mark_groups=True, palette=ig.RainbowPalette(n=3), edge_width=0.5, target=ax)
    plt.axis('off')
    plt.show()
    
