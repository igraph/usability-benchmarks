import igraph as ig
import matplotlib.pyplot as plt

g = ig.Graph.Famous("Zachary")

# Use edge betweenness to detect communities
communities = g.community_edge_betweenness()

# Convert the VertexDendrogram into a VertexCover for plotting
communities = communities.as_clustering().as_cover()
num_communities = len(communities)

# Visualize the communities simply with vertex cover coloring
fig, ax = plt.subplots()
ig.plot(
    communities,
    mark_groups=True,
    palette=ig.RainbowPalette(n=num_communities),
    edge_width=0.5,
    target=ax,
)
ax.legend(range(num_communities))
plt.show()
