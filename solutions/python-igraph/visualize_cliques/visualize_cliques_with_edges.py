import igraph as ig
import matplotlib.pyplot as plt

g = ig.Graph.Famous('Zachary')
g.vs['label'] = [i + 1 for i in range(g.vcount())]
g.vs['size'] = 0
g.es['width'] = 0.3
g.vs['color'] = 'yellow'

fig, ax = plt.subplots(3, 4)
cliques = g.cliques(4, 4)
for i in range(len(cliques)):
    g.vs[cliques[i]]['size'] = 0.5
    for j in range(len(cliques[0]) - 1):
        for k in range(j + 1, len(cliques[0])):
            g.es[g.get_eid(cliques[i][j], cliques[i][k])]['color'] = 'red'
            g.es[g.get_eid(cliques[i][j], cliques[i][k])]['width'] = 1
            
    ig.plot(ig.VertexCover(g, [cliques[i]]), mark_groups=True, target=ax[i // 4, i % 4])
    
    g.vs['size'] = 0
    g.es['color'] = 'black'
    g.es['width'] = 0.3
    
plt.axis('off')
plt.show()

