(* ::Package:: *)

(* ::Section:: *)
(*Vertex attributes*)


(* ::Text:: *)
(*Load package and create example graph.*)


Needs["IGraphM`"]
g = IGFamousGraph["Zachary"];


(* ::Text:: *)
(*Set the Distance vertex attribute to be the distance from vertex 1.*)


g // IGVertexMap[#&, "Distance" -> (First@IGDistanceMatrix[#, {1}]&)]


(* ::Section:: *)
(*Connected components*)


(* ::Text:: *)
(*Load package and create example graph.*)


Needs["IGraphM`"]
g = IGErdosRenyiGameGNM[50, 30];


(* ::Text:: *)
(*This solution does not make any use of IGraph/M. It is fully based on Mathematica built-ins.*)


connectedComponents[g_?GraphQ] :=
	Module[{f},
		Set @@@ Map[f] /@ Join[EdgeList[g], Reverse /@ EdgeList[g]];
		GatherBy[VertexList[g], f]
	]


connectedComponents[g]


(* ::Section:: *)
(*Number of edges incident to each vertex*)


(* ::Text:: *)
(*Create example graph.*)


n = 6; m = 12;
g = Graph[Range[n], UndirectedEdge @@@ RandomInteger[{1, n}, {m, 2}], VertexLabels->Automatic]


(* ::Text:: *)
(*This solution does not use IGraph/M. It fully relies on Mathematica built-ins. This is an efficient solutions because AdjacencyMatrix returns a sparse matrix.*)


Total@AdjacencyMatrix[g]


(* ::Section:: *)
(*Mean neighbour degree*)


(* ::Text:: *)
(*Load package and create example graph.*)


Needs["IGraphM`"]
g = IGFamousGraph["Zachary"];


meanNeighborDegree[g_?GraphQ, k_] := 
	With[{deg = IGVertexAssociate[VertexDegree][g]},
		Mean@Lookup[deg, AdjacencyList[g, #, k]]& /@ VertexList[g]
	]


meanNeighborDegree[g, 1]
meanNeighborDegree[g, 2]


(* ::Section:: *)
(*Edge multiplicities*)


(* ::Text:: *)
(*Load package and create example graph.*)


Needs["IGraphM`"]

n = 6; m = 12;
g = Graph[Range[n], UndirectedEdge @@@ RandomInteger[{1, n}, {m, 2}], VertexLabels->Automatic]


IGWeightedSimpleGraph[g, SelfLoops->True]


(* ::Section:: *)
(*Identifying multi-edges*)


n = 20; m = 60;
g = Graph[Range[n], UndirectedEdge @@@ RandomInteger[{1, n}, {m, 2}], VertexLabels->Automatic]


(* ::Text:: *)
(*Get vertex pairs between which there is more than one edge. This solution is written to work for edge-tagged graphs.*)


Cases[
	Gather[Sort /@ EdgeList[g][[All, {1,2}]]],
	{e_, __} :> e
]


(* ::Text:: *)
(*Groups of edge indices for edges between the same vertex pairs:*)


Values@GroupBy[
	Transpose@List[
		Range@EdgeCount[g],
		Sort /@ EdgeList[g][[All, {1,2}]]
	],
	Last -> First
]


(* ::Text:: *)
(*Group edges using tagged edge representation:*)


tg = EdgeTaggedGraph[g];
GatherBy[
	EdgeList[tg],
	Sort[#[[1;;2]]]&
]


(* ::Text:: *)
(*Merge parallel edges and count multiplicities. Preserve self-loops.*)


IGWeightedSimpleGraph[g, SelfLoops->True]
