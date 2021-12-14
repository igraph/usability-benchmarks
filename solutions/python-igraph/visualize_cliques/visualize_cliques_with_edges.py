import igraph as ig
import matplotlib.pyplot as plt

# Set up the base Zachary Karate Club graph
g = ig.Graph.Famous("Zachary")
g.vs["label"] = [i + 1 for i in range(g.vcount())]
base_visual_style = {
    "vertex_size": 0,
    "vertex_color": "white",
    "edge_width": 0.3
}

# Generate a sorted list of all eleven 4-cliques
cliques = g.cliques(4, 4)

# Generate 11 subplots, each highlighting one 4-clique and its edges
clique_visual_styles = []
for i in range(len(cliques)):
    clique_visual_styles.append(dict(base_visual_style))
    vertex_colors = ["white"] * g.vcount()
    edge_colors = ["black"] * g.ecount()
    vertex_sizes = [0] * g.vcount()
    edge_widths = [0.3] * g.ecount()

    for j in cliques[i]:
        vertex_colors[j] = "yellow"
        vertex_sizes[j] = 0.6
        
    for j in range(len(cliques[0]) - 1):
        for k in range(j + 1, len(cliques[0])):
            edge_colors[g.get_eid(cliques[i][j], cliques[i][k])] = "red"
            edge_widths[g.get_eid(cliques[i][j], cliques[i][k])] = 1

    clique_visual_styles[i]["vertex_color"] = list(vertex_colors)
    clique_visual_styles[i]["vertex_size"] = list(vertex_sizes)
    clique_visual_styles[i]["edge_color"] = list(edge_colors)
    clique_visual_styles[i]["edge_width"] = list(edge_widths)

# Create 3 rows of 4 subplots
fig, ax = plt.subplots(3, 4)
for i in range(len(clique_visual_styles)):
    ig.plot(ig.VertexCover(g, [cliques[i]]), mark_groups=True, target=ax[i // 4, i % 4], **clique_visual_styles[i])
plt.show()

