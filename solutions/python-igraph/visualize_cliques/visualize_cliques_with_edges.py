import igraph as ig
import matplotlib.pyplot as plt

g = ig.Graph.Famous('Zachary')

fig, ax = plt.subplots(3, 4)
cliques = g.cliques(4, 4)
for i, clique in enumerate(cliques):
    g.vs['color'] = 'yellow'
    g.es['color'] = 'black'
    g.es['width'] = 0.3
    
    g.vs[clique]['color'] = 'red'
    clique_edges = g.es.select(_within=clique)
    clique_edges['color'] = 'red'
    clique_edges['width'] = 1
            
    ig.plot(ig.VertexCover(g, [clique]), mark_groups=True, palette=ig.RainbowPalette(), target=ax[i // 4, i % 4])
  
plt.axis('off')
plt.show()
