import igraph as ig

def weakly_connected_components(g : ig.Graph):
    parent = sorted(list(range(len(g.vs))))

    def find(x):
        if x != parent[x]:
            parent[x] = find(parent[x])
        return parent[x]
    
    def union(x, y):
        parentx, parenty = find(x), find(y)
        if parentx != parenty:
            parent[parentx] = parenty

    for e in g.es:
        union(e.source, e.target)
    connected_components = dict()
    for i in range(20):
        p = find(parent[i])
        if p not in connected_components:
            connected_components[p] = []
        connected_components[p].append(i)
    return connected_components

g = ig.Graph.Erdos_Renyi(20, m=8, directed=False, loops=False)
wcc = weakly_connected_components(g)
colors = ["aliceblue", "antiquewhite", "aqua", "aquamarine", "beige", "navyblue", "orchid", "blue", "blueviolet", "brown", "cadetblue", "chartreuse", "cornflowerblue", "gold", "darkolivegreen", "darksalmon", "firebrick", "goldenrod", "khaki"]
for a in range(len(wcc.keys())):
    i = list(wcc.keys())[a]
    for j in wcc[i]:
        g.vs[j]["color"] = colors[a]
ig.plot(g, target='./solutions/python-igraph/connected_components/graph.png')