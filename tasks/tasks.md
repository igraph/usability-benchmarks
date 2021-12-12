
Since not all libraries have all functionality necessary for each task, prerequisites are noted whenever they may not be obvious.

# Basic tasks and algorithms

## Vertex attributes

Add a new numerical vertex attribute which is the distance of the respective vertex from a given vertex `v`.


## Connected components

Implement finding (weakly) connected components from scratch. Do not use the library's built-in connected component finder.


## Number of edges incident to each vertex

Write a function that computes how many edges are incident to each vertex. In undirected graphs with self-loops, this is not the same as the degree. An isolated vertex with a self-loop has one incident edge, but it has degree 2.

Create a random graph with `n` vertices and `m` edges, allowing self-loops. Demonstrate the function on it. `n=6`, `m=12` produces a small example that is easy to check visually.


## Mean neighbour degree

 1. Write a function to compute the mean degree of the `k`th order neighbourhood of each vertex, excluding the degree of the vertex itself. Write this function from scratch.
 2. Compute the result of the the karate club network for `k=1` and `k=2`.


## Edge multiplicities

 1. Create a random multigraph, allowing self-loops as well. Ensure it has some parallel self-loops.
 2. Create a new graph where parallel edges are merged. Self-loops should not be removed. The edge multiplicities of the original graph should be stored in an edge attribute of the new graph.


## Identifying multi-edges

 1. Create a random multigraph with `n` vertices and `m` edges. Both multi-edges and self-loops are allowed.
 2. Identify vertex pairs between which there is more than one edge.
 3. Create a list of edge-groups corresponding to the above. In igraph, edges should be identified based on their index, thus the output is a list of edge index lists.
 4. Simplify the multi-edges, but keep self-loops. The resulting graph should have edge weights which are equal to the edge multiplicity in the original graph.

Suggested parameters: `n=20`, `m=60`.


## Block-cut tree

**Prerequisite:** Functionality for biconnected components.

Compute the block-cut tree of a graph, and visualize it. The tree should be stored as a graph object. It should have a vertex attribute which stores the cut vertex names, or the blocks (as vertex sets). When visualizing the tree, use this attribute for labelling the vertices.


## Local complement

Create a new graph `g2` based on an existing one `g` so that two neighbours of a vertex `v` in `g` will be connected in `g2` precisely when they are not connected in `g`.


## Orient a tree

Write a function that:

 - Takes an undirected tree as input, as well as a vertex. This vertex will be the root.
 - Outputs a directed tree where all edges are oriented away from the root.


## Strahler stream order

Write a function that computes [the Strahler stream order](https://en.wikipedia.org/wiki/Strahler_number) of each node in a directed out-tree. The function should verify that the input is valid (i.e. it is a directed tree where all edges point away from the root).


# Network analysis

## Clique percolation

**Prerequisite:** Clique finder.

The clique percolation method can identify overlapping communities. It proceeds like this:

 - Find all `k`-cliques.
 - Two cliques that overlap in at least `k-1` vertices are considered connected.
 - A community is a connected component of cliques according to the above connectivity rule.

The task:

 1. Write a function that implements clique percolation for arbitrary `k`.
 2. Compute the result for `k=3` for the karate club network.
 3. Visualize it, making sure to handle overlapping communities correctly.


## Assortativity in connected graphs

 1. Generate `k` random graphs with `n` vertices and `m` edges each, and record their assortativities in a list.
 2. The same as above, but this time only generate connected graphs. Use rejection sampling.
 3. Compare the mean and standard deviation of the assortativity values in these two cases. Optionally, plot and compare histograms.

Suggested parameters: `k=10000, n=30, m=33`.


## Giant component

Create a plot illustrating the giant component phase transition in Erdős-Rényi graphs on `n` vertices. The horizontal axis must be the (edge count) / (vertex count) ratio and range from 0 to 2. The vertical axis must be the fraction of vertices contained within the giant component.

Suggested parameter: `n=10000`


## Targeted attack

Create an undirected preferential attachment graph on `n` vertices, adding `k` edges in each step. Keep removing the highest betweenness vertex until only `n/2` vertices remain. At each step, record the fraction of vertices in the largest component. Plot this value as a function of the number of steps taken.

Suggested parameters: `n=200`, `k=3`.


## Voronoi partitions

 1. Create a _connected_ geometric random graph with `n=100` vertices, cutoff `0.13`. Use rejection sampling to ensure connectedness.
 2. Weight the edges by their geometric length.
 3. Select three vertices at random. These are our centrepoints.
 4. For each centrepoint, identify the group of vertices which are closer to it (in graph distance) than to the others. This gives a partitioning of the vertices. Plot the partitions. The plot should use the original coordinates of the geometric random graph.


## Betweenness spanning tree

Highlight (visualize) the minimum and maximum spanning trees in the karate club network using edge betweenness as the edge weights.


## Clustering tree

**Prerequisite:** Any hierarchical clustering method.

Construct a rooted directed tree representing the hierarchical clustering obtained from a community detection method of your choice. The output should be a standard graph object.


# Graph theory

## Vertex cover

**Prerequisite:** Clique finder.

Find a single minimum vertex cover of a graph based on cliques. Do not use any built-in vertex cover function.

Implementation hints: Every minimum vertex cover is a complement of a maximum independent vertex set. An independent vertex set is a clique in the complement graph.


## Isomorphic duplicates

Generate 1000 random graphs with 5 vertices and 5 edges. Remove isomorphic duplicates from the result list.


# Visualization

## Community visualization

**Prerequisite:** Any community finding algorithm.

 1. Identify communities in the Zachary karate club network.
 2. Visualize the communities in whichever way is simplest.
 3. Add a legend to the plot, showing the index of each community (e.g. a numeric label for each community colour).
 4. Visualize the communities specifically by assigning a different colour to each community (both to edges and vertices). Edges running between communities should be grey.


## Visualize multigraphs

 1. Generate a small graph which has multi-edges and self-loops, including multiple self-loops on some of the vertices. It should have no more than 26 edges in total.
 2. Store the letters of the alphabet, A, B, C, etc., in an edge attribute named `letter`.
 3. Assign arbitrary numeric weights to each edge within the interval [0, 1].
 4. Visualize the graph such that the edges would be coloured according to their weight and a gradient colourmap of your choice. The edges should be labelled using their letter.

This tasks tests that:

 - Multigraphs can be plotted correctly, including multiple self-loops on a vertex.
 - Edge labels are positioned correctly in multigraphs (i.e. they are on a curved edge, not on the straight line connecting their endpoints).
 - The visualization function distinguishes between parallel edges correctly, i.e. each edge has its correct label and colour (not the label or colour of its parallel neighbour).


## Visualize cliques

**Prerequisite:** Clique finder.

 1. Identify all 4-cliques in the Zachary karate club network. There are 11 in total.
 2. Create a figure with 11 subplots, each showing one clique highlighted in the graph in whatever way.
 3. Same as above, but this time the highlighting should be done by using a separate style for both the vertices and edges contained within the clique.


## Visualize vertex properties using size

Visualize the eigenvector centrality in the Zachary karate club network by making the area (not radius) of each vertex be proportional to its centrality.


## Visualize vertex and edge properties using colour

 1. Generate a random geometric graph with 100 vertices and cutoff 0.15.
 2. Set the edge weights to be proportional to the geometric distance of their endpoints.
 3. Visualize the graph using the actual vertex positions.
 4. Colour vertices and edges based on their betweenness value, using a gradient scale. Vertex and edge betweenness should be scaled to their respective colour scales separately.
 5. Same visualization as above, but now vertex and edge betweenness values should be scaled together (i.e. they should be comparable to each other). Add a colour bar as plot legend. Label the plot with the description "Edge and vertex betweenness".


## Custom edge shapes

This is to test the flexibility of the visualization API.

 1. Plot a small graph with all edges being drawn as S-shaped squiggly lines.
 2. Plot a 3-path so that the first edge is squiggly and the second is straight. The ability to have different shapes for different edges is important for some layouts, such as Sugiyama.
 3. Plot a directed graph so that reciprocal edges are curved (i.e. they do not overlap) but unidirectional edges are straight.
 4. Plot a directed graph so that reciprocal edges are represented as straight bidirectional arrows, while unidirectional edges are represented as straight unidirectional arrows.
 5. Plot this multigraph so that simple edges are straight and multi-edges are curved: `1-2, 1-2, 1-2, 2-3, 3-1, 3-1`.
 6. Plot a bipartite graph in the style of a Sankey diagram. The two partitions should be laid out in two columns, left and right. The edges should be thick smoothly curving lines, with a horizontal tangent at both ends.
 7. Visualize a graph without drawing the edges. Vertices should be drawn as usual, retaining all possibilities for vertex shape customization.


## Custom vertex shapes

This is to test the flexibility of the visualization API.

 1. Plot a directed 3-cycle graph so that the vertices are a disk, a square and a triangle.
 2. Plot a two-path with the first vertex being a tall rectangle and the second a wide rectangle. Both vertices should be labelled underneath. Is label positioning taking into account the vertex shape?
 3. Plot a graph without drawing the vertices. Only the edges should be drawn. The edges should touch at their endpoints.
 4. Plot a simple food web so that the vertices are images of the plants and animals. Choose your food web so that it is acyclic, and use the Sugiyama layout. Example: https://www.azolifesciences.com/article/Threats-in-the-Food-Chain.aspx


## Implement hierarchical edge bundling

**Prerequisite:** Any method for hierarchical clustering of vertices.

The basic idea of hierarchical edge bundling is:

 - Compute a hierarchical clustering of the graph
 - Lay out the vertices on a circle, and use a radial layout of the clustering tree to guide edge shapes

There are many descriptions of the technique on the internet. @szhorvat provides a description as well as a Mathematica implementation here:  https://mathematica.stackexchange.com/q/55367/12


## Visualize an Eulerian path

Take an Eulerian graph (for example, octaherdal graph), and visualize an Eulerian cycle by colouring the edges along the cycle using a gradient colour map.


## Graph atlas, connected components

This is inspired by a [networkx example](https://networkx.org/documentation/stable/auto_examples/graphviz_layout/plot_atlas.html#sphx-glr-auto-examples-graphviz-layout-plot-atlas-py)

 1. Take graphs 53 through 208 from the graph atlas. These are all 6-vertex graphs.
 2. Filter the connected ones.
 3. Create their disjoint union, obtaining a many-component graph.
 4. Plot the graph so that each connected component has a different colour.


## Highlight path, dashed lines

Find a longest shortest path (diameter path) in the karate club network. Visualize this graph so that this path is shown:

 - thicker than other edges
 - in a different colour
 - drawn as a dashed line


## Biconnected components

Visualize biconnected components similarly to this example:

https://en.wikipedia.org/wiki/Biconnected_component#/media/File:Graph-Biconnected-Components.svg

 - Each component will have its assigned colour
 - Edges within a component are drawn using that colour
 - Vertices are drawn using all the colours of the components they belong to.
