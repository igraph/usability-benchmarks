## Tasks ref: https://github.com/igraph/usability-benchmarks/blob/master/tasks/tasks.md

                                        # Tasks headers are marked with ### ... ###
                                        # Task text marked with ## ...
                                        # My comments marked with # ...

library(igraph)

### Vertex Attributes ###

## Add a new numerical vertex attribute which is the distance of the respective vertex from a given vertex v.
                                        # Generate an arbitrary random graph, in this case with
                                        # 20 nodes and a density of 0.25
g <- sample_gnp(n = 20, p = .25)
                                        # Choose an arbitrary node
selectnode <- 1
                                        # The igraph function distances takes the graph, `g`, and the
                                        # chosen vertex, `selectnode`, and outputs the distance to
                                        # all other nodes `V(g)`. That output can be stored directly
                                        # as a vertex attribute.
V(g)$d <- distances(g, v = V(g)[selectnode], to = V(g))
                                        # To show that it worked:
plot(g, vertex.color = V(g)$d)
V(g)$d

### Connected Components ###

## Implement finding (weakly) connected components from scratch. Do not use the library's built-in connected component finder.

                                        # Solution based on depth first search algorithm in Cormen et
                                        # al 2009 Intro to Algorithms 3rd Ed

                                        # Generate a simple network that is easy to visualize.
n1 <- 5
n2 <- 5
g <- sample_gnm(n1, n1, TRUE) %>% set_vertex_attr("label", value = 1:n1)
h <- sample_gnm(n2, n2, TRUE) %>% set_vertex_attr("label", value = (n1+1):(n1+n2))
G <- g + h
plot(G)

                                        # This is an implementation of depth first search.
                                        # igraph has its own implementation as well, called `dfs`,
                                        # but because the task asked for a solution "from scratch",
                                        # here's a way to implement DFS in igraph. DFS is useful here
                                        # because it "discovers" the graph. The plan is to make the
                                        # graph collection, G, undirected so that DFS "discovers"
                                        # weakly connected components.

                                        # There is one change from the textbook algorithm to support
                                        # this particular task, marked below. 
dfs_visit <- function(g, u) {
    require(igraph)
                                        # The recursive part
    g$time <- g$time + 1
    V(g)$d[u] <- g$time
    V(g)$color[u] <- "gray"
    nbs <- neighbors(g, u, "out")
    for(v in nbs) {
        if(V(g)$color[v] == "white") {
            V(g)$p[v] <- as.numeric(u)
            g <- dfs_visit(g, V(g)[v])
        }
    }
                                        # Instead of always coloring a node black, use the chosen
                                        # mark color
    V(g)$color[u] <- g$mark_color
    g$time <- g$time + 1
    V(g)$f[u] <- g$time
    return(g)
}

depth_first_search <- function(g) {
    require(igraph)
    V(g)$color <- "white"
    V(g)$p <- numeric(length(V(g))) # predecessor
    V(g)$d <- numeric(length(V(g))) # discovery time
    V(g)$f <- numeric(length(V(g))) # finished time
    g$time <- 0 # time
                                        # Change from textbook example:
                                        # change mark color if the algorithm finishes DFS from a node.
                                        # This mark color can then be used to identify connected
                                        # components once the graph has been forced to be undirected.
                                        # This procedure won't always work as expected in a directed
                                        # graph, for example if the graph is a line, but that
                                        # limitation isn't critical here.
    g$mark_color <- 1 
    for(i in 1:vcount(g)) {
        if(V(g)$color[i] == "white") {
            g <- dfs_visit(g, V(g)[i])
        }
                                        # Increment mark color
        g$mark_color <- g$mark_color + 1
    }
    return(g)
}

weakly_connected_components <- function(G) {
    require(igraph)
                                        # Connected components in an undirected graph correspond to
                                        # weakly connected components in a directed graph.
    g <- as.undirected(G, mode = "collapse")
    g <- depth_first_search(g)
                                        # Now the mark color tells us the components.
    comps <- unique(V(g)$color)
    components <- list()
    for(i in 1:length(comps)) components[[i]] <- V(g)[which(V(g)$color == comps[i])]
                                        # Return the components, not the graph; not the original
                                        # graph remains unaltered.
    return(components)
}
    
weakly_connected_components(G)

### Number of edges incident to each vertex ###

## Write a function that computes how many edges are incident to each vertex. In undirected graphs with self-loops, this is not the same as the degree. An isolated vertex with a self-loop has one incident edge, but it has degree 2.
count_incident_edges <- function(g) {
    require(igraph)
                                        # Get all edges in a matrix format
    es <- ends(g, es = E(g))
                                        # And the vertex labels
    vs <- V(g)
                                        # Count the number of edges connected to a node
                                        # This function works by checking, for each vertex,
                                        # how many rows in `es` have at least one instance of vertex
                                        # label `v`.
    edge_count <- sapply(vs, function(v) sum(apply(es, 1, function(e) v %in% e)))
    names(edge_count) <- vs
    return(edge_count)
}

## Create a random graph with n vertices and m edges, allowing self-loops. Demonstrate the function on it. n=6, m=12 produces a small example that is easy to check visually.
g <- sample_gnm(n = 6, m = 12, directed = FALSE, loops = TRUE)
plot(g)
                                        # This is a numeric vector, the names of which are the vertex
                                        # labels and the values are the counts of the edges incident
                                        # to the vertex.
count_incident_edges(g)

### Mean neighbour degree ###

## Write a function to compute the mean degree of the kth order neighbourhood of each vertex, excluding the degree of the vertex itself. Write this function from scratch.
mean_neighbor_degree <- function(g, v, k) {
    require(igraph)
                                        # The neighborhood function allows you to specify an order,
                                        # which neighbors does not, but includes the node `v`
    nbs <- neighborhood(g, order = k, nodes = v, mode = "all")[[1]]
                                        # so we can just take it out.
    nbs <- nbs[which(nbs != v)]
                                        # Once we have the order of neighbors we want, we can collect
                                        # the degree for those nodes
    nbs_k <- sapply(nbs, function(n) degree(g, n))
                                        # and return the mean
    return(mean(nbs_k))
}

                                        # It is easier to make the above calculations for every
                                        # node if we define another sapply() call.
all_mean_neighbor_degree <- function(g, k) {
    sapply(V(g), function(v) mean_neighbor_degree(g, v, k))
}

                                        # The Karate graph is also in the package "igraphdata",
                                        # which has other useful/historical graphs, and in
                                        # "networkdata", by David Schoch and also including the
                                        # weighted version, available here:
                                        # http://networkdata.schochastics.net/
g <- make_graph("Zachary")
                                        # Choose an arbitrary node
v <- V(g)[7]
                                        # Mean neighbor degree for nodes one step away
mean_neighbor_degree(g, v = v, k = 1)
                                        # and two steps away.
mean_neighbor_degree(g, v = v, k = 2)

## Compute the result of the the karate club network for k=1 and k=2.
all_mean_neighbor_degree(g, k = 1)
all_mean_neighbor_degree(g, k = 2)


### Edge multiplicities ###

## Create a random multigraph, allowing self-loops as well. Ensure it has some parallel self-loops.
random_multigraph <- function(n, m, directed = TRUE) {#, req.parallel = TRUE, max.iter = 10) {
    "Create a multigraph with n nodes and m edges."
                                        # The basic sample_(gnm | gnp) don't want to create
                                        # multigraphs, so this function will make a basic GNM
                                        # multigraph with a specified number of nodes and edges.

                                        # A vector for nodes
    nodes <- seq(1, n)
                                        # and a storage list for edges.
    edges <- vector("list", m)
                                        # For each expected edge, 
    for(i in 1:m) {
                                        # choose a node at random 
        u <- sample(nodes, 1)
                                        # and another 
        v <- sample(nodes, 1)
                                        # make an edge
        edge <- c(u, v)
                                        # and add that edge to our storage list.
        edges[[i]] <- edge
    }
                                        # Collapse the list into a vector
    edges <- do.call(c, edges)
                                        # and make the graph.
    g <- graph(edges)
                                        # Make undirected if desired.
                                        # The argument mode = "each" is necessary to preserve the
                                        # multigraph.
    if(!directed) g <- as.undirected(g, mode = "each")
                                        # And return the graph.
    return(g)
}

sample_random_multigraph <- function(n, m,
                                     directed = TRUE,
                                     req.parallel = TRUE,
                                     req.self.parallel = TRUE,
                                     max.iter = 10) {
    require(igraph)
                                        # This function calls `random_multigraph` until it produces a
                                        # graph meeting the requirements specified in the arguments.
    if(req.parallel | req.self.parallel) {
        test <- FALSE
                                        # Check for requirements in a `for` loop to prevent infinite
                                        # loops.
        for(i in 1:max.iter) {
                                        # Make a trial graph
            g <- random_multigraph(n, m, directed = directed)
                                        # and the edgelist that goes a long with it, for checking
            el <- as_edgelist(g)
                                        # Parallel edges will show up as duplicated rows in `el`
                                        # If there are no duplicated rows, try again.
            if(!any(duplicated(el))) {
                next
            } else if(!req.self.parallel) {
                                        # If self parallel edges are not required, return the graph
                return(g)
            }
            
                                        # Self parallel edges are a special case.
            if(req.self.parallel) {
                                        # This returns all copies of duplicated rows
                dupls <- el[duplicated(el) | duplicated(el, fromLast = TRUE), ]
                                        # So the test for self-parallel is whether or not there are
                                        # any rows with two entries of the same node label.
                test <- any(dupls[, 1] == dupls[, 2])
                                        # If there are double entries the requirements are satisfied.
                if(test) break
            }
        }
                                        # If no graph could be made satifying the requirements, tell
                                        # the user
        if(!test) return("Could not make required graph in max iterations.")
    } else {
                                        # If additional specifications aren't needed, just make the
                                        # graph
        g <- random_multigraph(n, m, directed = directed)
    }
                                        # If conditions are satisfied, return the graph
    return(g)

}

g <- sample_random_multigraph(5, 10, directed = FALSE)
    
## Create a new graph where parallel edges are merged. Self-loops should not be removed. The edge multiplicities of the original graph should be stored in an edge attribute of the new graph.
simplify_mod <- function(g) {
                                        # There is a function that simplifies graphs, but it needs
                                        # support for this task: an edge attribute that counts edges
    E(g)$edge_count <- 1
    g <- simplify(g, remove.multiple = TRUE, remove.loops = FALSE, edge.attr.comb = "sum")
    return(g)
}

h <- simplify_mod(g)

E(g)
E(h)
E(h)$edge_count

### Identifying Multi-Edges ###

## 1. Create a random multigraph with n vertices and m edges. Both multi-edges and self-loops are allowed.
                                        # Using the function above...
## Suggested parameters: n=20, m=60.
g <- random_multigraph(n = 20, m = 60, directed = FALSE)

## 2. Identify vertex pairs between which there is more than one edge.
                                        # Using the same strategy with duplicated rows in an array
el <- as_edgelist(g)
dupls <- el[duplicated(el), ]
                                        # but putting duplicate node pairs in a list.
dupls <- lapply(1:nrow(dupls), function(i) dupls[i, ])

## 3. Create a list of edge-groups corresponding to the above. In igraph, edges should be identified based on their index, thus the output is a list of edge index lists.
                                        # Applying a function to the list of node pairs,
edge_groups <- lapply(dupls, function(pair) {
                                        # Get all edges for each node separately
    edges <- incident_edges(graph = g, v = pair)
                                        # Multi edges between those nodes show up in both values
                                        # of the incident_edges() output list.
    multi <- edges[[1]][which(edges[[1]] %in% edges[[2]])]
                                        # Return those edges.
    multi
})

## 4. Simplify the multi-edges, but keep self-loops. The resulting graph should have edge weights which are equal to the edge multiplicity in the original graph.
                                        # We can use the function written above
h <- simplify_mod(g)
                                        # and check to make sure it's correct.
edge_groups
E(h)[which(E(h)$edge_count > 1)]
E(h)[which(E(h)$edge_count > 1)]$edge_count


### Block-cut Tree ###

                                        # Generate a random graph
g <- sample_gnm(10, 14)
                                        # and count the number of biconnected components.
biconnected_components(g)$no

## Compute the block-cut tree of a graph, and visualize it. The tree should be stored as a graph object. It should have a vertex attribute which stores the cut vertex names, or the blocks (as vertex sets). When visualizing the tree, use this attribute for labelling the vertices.
                                        # The biconnected_components() function has the information
                                        # we need.
bic <- biconnected_components(g)
                                        # The $components list has the vertices that go in each
                                        # component
comps <- bic$components
                                        # Label for convenience
names(comps) <- LETTERS[1:length(comps)]
                                        # $articulation_points has which nodes belong in multiple
                                        # components
cut_v <- bic$articulation_points

                                        # Label for convenience
V(g)$color <- ifelse(V(g) %in% cut_v, "dodgerblue", "orange")
V(g)$comp <- sapply(V(g), function(v) {
    which_comps <- names(comps)[which(sapply(comps, function(x) v %in% x))]
    paste(which_comps, collapse = ", ")
})

                                        # Make the block-cut tree as a graph
bic_tree <- unfold_tree(g, roots = cut_v)
h <- bic_tree$tree
                                        # and keep track of the original node labels
idx <- bic_tree$vertex_index
                                        # Label with the information gathered above
V(h)$color <- sapply(idx, function(v) V(g)[v]$color)
V(h)$comp <- sapply(idx, function(v) V(g)[v]$comp)

par(mfrow = c(2, 2))
coords1 <- layout_nicely(g)
plot(g, layout = coords1)
plot(g, vertex.label = V(g)$comp, layout = coords1)
coords2 <- layout_nicely(h)
plot(h, vertex.label = idx, layout = coords2)
plot(h, vertex.label = V(h)$comp, layout = coords2)

### Local complement ###

## Create a new graph g2 based on an existing one g so that two neighbours of a vertex v in g will be connected in g2 precisely when they are not connected in g.
                                        # Make a graph
g <- sample_gnm(5, 6)
                                        # in a circle, so connections are easy to see.
coords <- layout_in_circle(g)
                                        # There is already an igraph function for this task.
g2 <- complementer(g)

par(mfrow = c(1, 2))
plot(g, layout = coords)
plot(g2, layout = coords)

### Vertex name handling in disjoint union ###

## Create the following two graphs, with the given vertex names:
## a-b-c-d-b
## e-f-d-e
                                        # This is done most easily with the graph_from_literal
                                        # function
g1 <- graph_from_literal(a-b, b-c, c-d, d-b)
g2 <- graph_from_literal(e-f, f-d, d-e)

## Produce their disjoint union, ...
G <- disjoint_union(g1, g2)

## then add a single edge between vertex d from the first graph and vertex d from the second graph. This task checks if it is easy to find out which vertex in the disjoint union corresponds to which vertex in the original graphs.
idx <- which(V(G)$name == "d")
G2 <- add_edges(G, c(V(G)[idx[1]], V(G)[idx[2]]))

par(mfrow = c(1, 2))
coords <- layout_nicely(G)
plot(G)
plot(G2, layout = coords)
