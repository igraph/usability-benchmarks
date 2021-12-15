import igraph as ig
import matplotlib.pyplot as plt

g = ig.Graph.Famous('Zachary')
cliques = g.cliques(4, 4)
fig, ax = plt.subplots(3, 4)
for i, clique in enumerate(cliques):
    ig.plot(ig.VertexCover(g, [clique]), mark_groups=True, palette=ig.RainbowPalette(), edge_width=0.5, target=ax[i // 4, i % 4])
plt.axis('off')
plt.show()
