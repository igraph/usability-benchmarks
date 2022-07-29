import igraph as ig
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

g = ig.Graph.Famous("Zachary")

# Use edge betweenness to detect communities
communities = g.community_edge_betweenness()

# Convert the VertexDendrogram into a VertexCover for plotting
communities = communities.as_clustering().as_cover()
num_communities = len(communities)

# Create a custom color legend
palette = ig.RainbowPalette(n=num_communities)
legend_elements = [
    Patch(facecolor=palette.get(i), edgecolor="k", label=i)
    for i in range(num_communities)
]

# Visualize the communities simply with vertex cover coloring
fig, ax = plt.subplots()
ig.plot(
    communities,
    mark_groups=True,
    palette=ig.RainbowPalette(n=num_communities),
    edge_width=0.5,
    target=ax,
)
ax.legend(handles=legend_elements)
plt.show()
