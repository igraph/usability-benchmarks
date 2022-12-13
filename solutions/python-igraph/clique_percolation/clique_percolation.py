import igraph as ig
import itertools
import matplotlib.pyplot as plt

EXAMPLE_K = 3

def clique_percolation(g, k, min_common_vertices):
    """
    Given a graph g and clique size k, return a list of communities found
    through clique percolation with at least min_common_vertices shared.
    A community is the connected component of cliques where two cliques overlap
    in at least min_common_vertices vertices.
    """
    if (min_common_vertices > k):
        raise ValueError(f"Error: min_common_vertices={min_common_vertices} must be greater than k={k}")

    # Cliques of size k in graph g
    k_cliques = g.cliques(k, k)

    # Each clique can be identified by it's index in 'cliques'
    cliques = [set(clique) for clique in k_cliques]
    num_cliques = len(cliques)

    # Get connected cliques
    # If the two cliques share at least min_common_vertices,
    # they are part of the same community
    cg_edges = [(c1, c2) for c1, c2 in itertools.combinations(range(num_cliques), 2) if len(cliques[c1].intersection(cliques[c2])) >= min_common_vertices]

    # Graph storing each clique as a vertex
    clique_graph = ig.Graph(cg_edges)

    # Merge cliques accordingly into communities
    cg_clusters = clique_graph.clusters()

    # Generate the list of clique communities as a list of set of vertices
    communities = []

    for g in cg_clusters:
        vertices = set()
        for c in g:
            vertices.update(cliques[c])
        communities.append(vertices)

    return communities

# Load in Zachary's karate club network
print("For the Zachary Karate Club graph\n")
g = ig.Graph.Famous("Zachary")

# Calculate cliques and communities
communities = clique_percolation(g, EXAMPLE_K, EXAMPLE_K-1)
print(f"Communities from cliques of size k={EXAMPLE_K}, with {EXAMPLE_K-1} shared vertices")
print(communities)

# Plot the graph with a vertex cover for each community
fig, ax = plt.subplots()
ig.plot(ig.VertexCover(g, communities), mark_groups=True, palette=ig.RainbowPalette(n=3), edge_width=0.5, target=ax)
plt.axis('off')
plt.show()

