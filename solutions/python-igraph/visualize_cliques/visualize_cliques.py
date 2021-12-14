import igraph as ig
import matplotlib.pyplot as plt

g = ig.Graph.Famous("Zachary")
cliques = g.cliques(4, 4)
fig, ax = plt.subplots(3, 4)
for i in range(len(cliques)):
    ig.plot(ig.VertexCover(g, [cliques[i]]), mark_groups=True, target=ax[i // 4, i % 4])
plt.show()

