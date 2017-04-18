# Minimum Spanning Trees

This repository holds the code for two approximation algorithms for partitioning graphs into balanced parts. It also contains code used to generate a specific set of planar graphs that were used to test these algorithms.

## Prerequisites
The spectural algorithm, Ncuts.cpp, requires the LAPACK library as well its c interface LAPACKE. Please refer to the LAPACKE documentation for instructions on how to properly set it up.

## Graph format

Graphs generated and used in this repo are undirected and use the following format:
```
p edge V E
e u1 v1  w1
e u2 v2  w2
e u3 v3  w3
...
```
Where V is the number of vertices in the graph, E the number of edges in the graph, u*x* the index of a vertex at one end of an edge, v*x* the index of a vertex at the other end of the edge and w*x* is the weight of the edge.
## Use
To compile the executables, simply run the following command:
```
make
```
To run the spectural clustering algorithm, run the following command replacing GRAPH_FILE and NUMBER_OF_BLOCKS with the name of the graphfile you wish to partition and the number of blocks you desire respectively.
```
./NCuts GRAPH_FILE  NUMBER_OF_BLOCKS
```
To run the dynamic programming algorithm, run the following command replacing GRAPH_FILE and NUMBER_OF_BLOCKS with the name of the graphfile you wish to partition and the number of blocks you desire respectively.
```
./DPCuts GRAPH_FILE  NUMBER_OF_BLOCKS
```

