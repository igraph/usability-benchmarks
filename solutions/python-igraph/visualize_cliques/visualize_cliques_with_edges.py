from igraph import Graph, plot
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # Set up the base Zachary Karate Club graph
    base = Graph.Famous("Zachary")
    base.vs["name"] = [i + 1 for i in range(base.vcount())]
    base.vs["label"] = base.vs["name"]
    base_layout = base.layout("kk")
    base_visual_style = {
        "vertex_size": 0,
        "vertex_color": "white",
        "edge_width": 0.3
    }

    # Generate a sorted list of all eleven 4-cliques
    cliques = base.cliques(4, 4)
    cliques.sort()
    
    # Generate 11 subplots, each highlighting one 4-clique and its edges
    clique_visual_styles = []
    for i in range(len(cliques)):
        clique_visual_styles.append(dict(base_visual_style))
        vertex_colors = ["white"] * base.vcount()
        edge_colors = ["black"] * base.ecount()
        vertex_sizes = [0] * base.vcount()
        edge_widths = [0.3] * base.ecount()

        for j in cliques[i]:
            vertex_colors[j] = "yellow"
            vertex_sizes[j] = 0.6
            
        for j in range(len(cliques[0]) - 1):
            for k in range(j + 1, len(cliques[0])):
                edge_colors[base.get_eid(cliques[i][j], cliques[i][k])] = "red"
                edge_widths[base.get_eid(cliques[i][j], cliques[i][k])] = 1

        clique_visual_styles[i]["vertex_color"] = list(vertex_colors)
        clique_visual_styles[i]["vertex_size"] = list(vertex_sizes)
        clique_visual_styles[i]["edge_color"] = list(edge_colors)
        clique_visual_styles[i]["edge_width"] = list(edge_widths)

    # Create 3 rows of 4 subplots
    fig, ax = plt.subplots(3, 4)
    plot(base, layout=base_layout, target=ax[0, 0], **base_visual_style)
    for i in range(len(clique_visual_styles)):
        plot(base, layout=base_layout, target=ax[(i + 1) // 4, (i + 1) % 4], **clique_visual_styles[i])
    plt.show()
