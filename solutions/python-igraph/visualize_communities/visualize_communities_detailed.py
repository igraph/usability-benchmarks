import igraph as ig
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

g = ig.Graph.Famous("Zachary")

# Use edge betweenness to detect communities
# and covert into a VertexClustering
communities = g.community_edge_betweenness()
communities = communities.as_clustering()
num_communities = len(communities)

# Color each vertex and edge based on its community membership
palette = ig.RainbowPalette(n=num_communities)
for i, community in enumerate(communities):
    g.vs[community]["color"] = i
    community_edges = g.es.select(_within=community)
    community_edges["color"] = i

# Create a custom color legend
legend_elements = [
    Patch(facecolor=palette.get(i), edgecolor="k", label=i)
    for i in range(num_communities)
]

# Plot with only vertex and edge coloring
fig, ax = plt.subplots()
ig.plot(communities, palette=palette, edge_width=1, target=ax)
ax.legend(handles=legend_elements)
plt.show()
