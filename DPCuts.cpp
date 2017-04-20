// C++ program for Kruskal's algorithm to find Minimum Spanning Tree
// of a given connected, undirected and weighted graph
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <utility>
#include "prettyPrint.h"

using namespace std;



struct Vertex
{
	int id;
	int mst_degree = 0;
};


// a structure to represent a weighted edge in graph
struct Edge
{
	static vector<Vertex> vertices;
	int src, dest, weight,x,y;
};


vector<Vertex> Edge :: vertices;

// a structure to represent a connected, undirected and weighted graph
struct Graph
{
	// V-> Number of vertices, E-> Number of edges
	int V, E;

	// graph is represented as an array of edges. Since the graph is
	// undirected, the edge from src to dest is also edge from dest
	// to src. Both are counted as 1 edge here.
	Edge* edge;
};


// Creates a graph with V vertices and E edges
Graph* createGraph(int V, int E)
{
	Graph* graph = (Graph*) malloc( sizeof(Graph) );
	graph->V = V;
	graph->E = E;

	for (int i = 0; i < V; ++i)
	{
		Vertex v;
		v.id = i;
		Edge::vertices.push_back(v);
	}
	graph->edge = (Edge*) malloc( graph->E * sizeof(Edge) );

	return graph;
}


// A structure to represent a subset for union-find
struct subset
{
	int parent;
	int rank;
};


// A utility function to find set of an element i
// (uses path compression technique)
int find(subset subsets[], int i)
{
	// find root and make root as parent of i (path compression)
	if (subsets[i].parent != i)
		subsets[i].parent = find(subsets, subsets[i].parent);

	return subsets[i].parent;
}


// A function that does union of two sets of x and y
// (uses union by rank)
void Union(subset subsets[], int x, int y)
{
	int xroot = find(subsets, x);
	int yroot = find(subsets, y);

	// Attach smaller rank tree under root of high rank tree
	// (Union by Rank)
	if (subsets[xroot].rank < subsets[yroot].rank)
		subsets[xroot].parent = yroot;
	else if (subsets[xroot].rank > subsets[yroot].rank)
		subsets[yroot].parent = xroot;

	// If ranks are same, then make one as root and increment
	// its rank by one
	else
	{
		subsets[yroot].parent = xroot;
		subsets[xroot].rank++;
	}
}


// Compare two edges according to their weights.
// Used in qsort() for sorting an array of edges
int myComp(const void* a, const void* b)
{
	Edge* a1 = (Edge*)a;
	Edge* b1 = (Edge*)b;
	if (a1->weight < b1->weight)
	{
		return -1;
	}
	else if (a1->weight > b1->weight)
	{
		return 1;
	}
	else
	{
		if (a1->src < b1->src)
		{
			return -1;
		}
		else if (a1->src > b1->src)
		{
			return 1;
		}
		else
		{
			if (a1->dest < b1->dest)
			{
				return -1;
			}
			else if (a1->dest > b1->dest)
			{
				return 1;
			}
			else
			{
				cout << "DUPLICATE EDGE; BAILING OUT" << endl;
				exit(0);
			}
		}
	}
}


void gimmeAdjListAndWeights(Graph *MST, vector<vector<int> > &adjList, vector<vector<int> > &distances)
{
	distances.resize(MST->V);
	for (int i = 0; i < distances.size(); ++i)
	{
		distances.at(i).resize(MST->V, -1);
		distances.at(i).at(i) = 0;
	}


	adjList.resize(MST->V);
	for (int i = 0; i < MST->V - 1; ++i)
	{
		// build up the adjacency list
		adjList.at(MST->edge[i].src).push_back(MST->edge[i].dest);
		adjList.at(MST->edge[i].dest).push_back(MST->edge[i].src);

		distances.at(MST->edge[i].src).at(MST->edge[i].dest) = MST->edge[i].weight;
		distances.at(MST->edge[i].dest).at(MST->edge[i].src) = MST->edge[i].weight;
	}

	// // TESTING:  print the adjacency list
	// cout << endl << endl << "adjList:  " << endl;
	// for (int i = 0; i < adjList.size(); ++i)
	// {
	// 	cout << i << ":  ";
	// 	for (int j = 0; j < adjList.at(i).size(); ++j)
	// 	{
	// 		cout << adjList.at(i).at(j) << "  ";
	// 	}
	// 	cout << endl;
	// }
	// cout << endl;
}

// TODO:  rewrite this to be iterative instead of recursive
void findDist(vector<vector<int> > &adjList, vector<vector<int> > &distances, int currentNode, vector<int> drapers, vector<int> &discovered)
{
	for (int i = 0; i < drapers.size(); ++i)
	{
		for (int j = 0; j < discovered.size(); ++j)
		{
			distances.at(drapers.at(i)).at(discovered.at(j)) = distances.at(currentNode).at(discovered.at(j)) + distances.at(currentNode).at(drapers.at(i));
			distances.at(discovered.at(j)).at(drapers.at(i)) = distances.at(currentNode).at(discovered.at(j)) + distances.at(currentNode).at(drapers.at(i));
		}
		discovered.push_back(drapers.at(i));


		vector<int> newDrapers; // vector of ints containing the nodes that "drape" from current draper
		// populate newDrapers:
		for (int k = 0; k < adjList.at(drapers.at(i)).size(); ++k)
		{
			if ( adjList.at(drapers.at(i)).at(k) != currentNode )
			{
				newDrapers.push_back( adjList.at(drapers.at(i)).at(k) );
			}
		}

		findDist(adjList, distances, drapers.at(i), newDrapers, discovered);
	}
}

void explore(vector<vector<int> > &adjList, vector<int> &path, vector<bool> &visited, int target)
{
	visited.at(path.back()) = true;
	if (path.back() == target)
	{
		return;
	}
	for (int i = 0; i < adjList.at(path.back()).size(); ++i)
	{
		if (!visited.at(adjList.at(path.back()).at(i)))
		{
			path.push_back(adjList.at(path.back()).at(i));
			explore(adjList, path, visited, target);
			if (path.back() == target)
			{
				return;
			}
			path.pop_back();
		}
	}
}

vector<int> findPath(vector<vector<int> > &adjList, vector<vector<int> > &distances)
{
	int max = -1, outer, inner;
	for (int i = 0; i < distances.size(); ++i)
	{
		for (int j = 0; j < distances.at(i).size(); ++j)
		{
			if (distances.at(i).at(j) > max)
			{
				max = distances.at(i).at(j);
				outer = i, inner = j;
			}
		}
	}

	// run dfs to find the path between furthest nodes
	vector<int> path(0);
	vector<bool> visited(adjList.size(), false);
	path.push_back(outer);
	explore(adjList, path, visited, inner);

	return path;
}

vector<int> miniFindPath(vector<vector<int> > &adjList, vector<vector<int> > &distances, bool* inBlock)
{
	int max = -1, outer, inner;
	for (int i = 0; i < distances.size(); ++i)
	{
		for (int j = 0; j < distances.at(i).size(); ++j)
		{
			if ( (distances.at(i).at(j) > max) && (inBlock[i]) && (inBlock[j]) )
			{
				max = distances.at(i).at(j);
				outer = i, inner = j;
			}
		}
	}

	// run dfs to find the path between furthest nodes
	vector<int> path(0);
	vector<bool> visited(adjList.size(), false);
	path.push_back(outer);
	explore(adjList, path, visited, inner);

	return path;
}

int branchSum(vector<vector<int> > &adjList, vector<vector<int> > &distances, int currentNode, int avoidNext, int avoidPrevious)
{
	int runningSum = 0;
	for (int i = 0; i < adjList.at(currentNode).size(); ++i)
	{
		if ( adjList.at(currentNode).at(i) != avoidNext && adjList.at(currentNode).at(i) != avoidPrevious )
		{
			runningSum += distances.at(currentNode).at(adjList.at(currentNode).at(i));
		}

		if (adjList.at( adjList.at(currentNode).at(i) ).size() > 1 && adjList.at(currentNode).at(i) != avoidNext && adjList.at(currentNode).at(i) != avoidPrevious )
		{
			runningSum += branchSum(adjList, distances, adjList.at(currentNode).at(i), currentNode, currentNode);
		}
	}
	return runningSum;
}

long long** partitionDPcuts(vector<vector<int> > &adjList, vector<vector<int> > &distances, vector<int> &path, int MSTweight, int numPartitions)
{
	// begin by creating the 2D array that will store the info about where to cut
	long long ** theFavor = new long long * [path.size()];
	for (int i = 0; i < path.size(); ++i)
	{
		// For a given node x after which we are considering a cut...
		// theFavor[x][0] = running best penalty at edge x
		// theFavor[x][1] = the previous node after which a cut should be made which results in the penalty at index 0 (best penalty achievable at this node)
		// theFavor[x][2] = total weight up to node x including branching off of node x but not including the weight between x and x + 1 (for calculating subtree weights)
		// theFavor[x][3] = weight of sub-portion of MST that occurs between nodes x and theFavor[x][1] (previous node that should be cut)
		// theFavor[x][4] = MSTweight - any edges (including edge x to x + 1) that have been cut to minimize penalty up to node x (for recalculating "ideal" block weight)
		theFavor[i] = new long long [5];
	}

	// NOTE:  the next line will segfault if your graph is too small.  The moral of the story:  use sufficiently big graphs
	int nextEdge = distances.at( path.at(0) ).at( path.at(1) ); // the weight of the edge on the path that comes after node cutAfterIndex
	int ideal = (MSTweight - nextEdge) / numPartitions; // target block weight
	int subMSTweight; // weight between previous edge we are considering cutting and current edge
	int penalty; // penalty incurred by making a given cut

	// set up the zero-index of theFavor
	theFavor[0][0] = ideal;
	theFavor[0][1] = 0;
	theFavor[0][2] = 0;
	theFavor[0][3] = 0;
	theFavor[0][4] = nextEdge;


	for (int i = 1; i < path.size(); ++i) // for each edge in the path
	{
		// assume for starters that the best way to cut at iteration i is to go straight from node 0 to node i
		if (i == path.size() - 1)
		{
			theFavor[i][2] = theFavor[i-1][2] + distances.at( path.at(i - 1) ).at( path.at(i) ); // note that there is necessarily no branching at the last node
			nextEdge = 0;
		}
		else
		{
			theFavor[i][2] = theFavor[i-1][2] + distances.at( path.at(i - 1) ).at( path.at(i) ) + branchSum(adjList, distances, path.at(i), path.at(i + 1), path.at(i - 1) );
			nextEdge = distances.at( path.at(i) ).at( path.at(i + 1) );
		}

		ideal = (MSTweight - nextEdge) / numPartitions; // figure out what the new ideal weight should be
		// set up theFavor with the info from the first iteration so we have something to compare the rest of the previous edges with
		theFavor[i][0] = abs(theFavor[i][2] - ideal); // penalty of makeing no cut
		theFavor[i][1] = i; // assume the cut should be made after node i (i.e. start by assuming the best cut is not to make a cut)...
		theFavor[i][3] = theFavor[i][2]; // ... in which case the weight of the subtree stored here should be the entire subtree "before" this node
		theFavor[i][4] = nextEdge; // ... and the only edge that we would cut in this case is the one after node i

		for (int j = 0; j < i; ++j) // for each edge that comes before the current edge (which is edge i)
		{
			// calculate penalty for this step
			ideal = (MSTweight - nextEdge - theFavor[j][4]) / numPartitions; // Update ideal
			subMSTweight = theFavor[i][2] - theFavor[j][2] - distances.at( path.at(j) ).at( path.at(j + 1) );
			penalty = theFavor[j][0] + abs(subMSTweight - ideal);
			if (penalty < 0)
			{
				cout << "overflow in penalty calculation!! " << endl;
				exit(0);
			}
			if ( penalty < theFavor[i][0] )
			{
				// reset the values of theFavor[i][0, 1, 3, and 4]
				theFavor[i][0] = penalty;
				theFavor[i][1] = j;
				theFavor[i][3] = subMSTweight;
				theFavor[i][4] = nextEdge + theFavor[j][4];
			}
		}
	}

	return theFavor;
}

// The main function to construct MST using Kruskal's algorithm
Graph *KruskalMST(Graph* graph, int numPartitions)
{
	int V = graph->V;
	// Edge result[V]; // This will store the resultant MST

	// let MST be an edge pointer holding the MST that we find
	Graph *MST = createGraph(graph->V, graph->V-1);

	int e = 0; // An index variable, used for result[]
	int i = 0; // An index variable, used for sorted edges

	// Step 1: Sort all the edges in non-decreasing order of their weight
	// If we are not allowed to change the given graph, we can create a copy of
	// array of edges
	qsort(graph->edge, graph->E, sizeof(graph->edge[0]), myComp);

	subset *subsets = (subset*) malloc( V * sizeof(subset) );

	// Create V subsets with single elements
	for (int v = 0; v < V; ++v)
	{
		subsets[v].parent = v;
		subsets[v].rank = 0;
	}

	vector<Edge> v;
	int prevEdgeWeight = -999999;
	// Number of edges to be taken is equal to V-1
	while (e < V - 1)
	{
		// Step 2: Pick the smallest edge. And increment the index
		// for next iteration
		Edge next_edge = graph->edge[i++];

		if(prevEdgeWeight != next_edge.weight && prevEdgeWeight > 0)
		{
			int minDegree;
			int currDegree;
			int index;
			Edge minEdge;
			minEdge.src = -1;
			while(v.size()!=0)
			{
				minDegree = 31231231;
				for (int j = 0; j < v.size(); ++j)
				{

					currDegree = Edge::vertices[v[j].src].mst_degree + Edge::vertices[v[j].dest].mst_degree; // maybe this should be max of src and dest and then minimize the max edge degree...
					if(currDegree < minDegree)
					{
						minDegree = currDegree;
						minEdge = v[j];
						index = j;
					}

				}
				v.erase(v.begin()+index);
				int x = find(subsets, minEdge.src);
				int y = find(subsets, minEdge.dest);
				minEdge.x = x;
				minEdge.y = y;
				// If including this edge does't cause cycle, include it
				// in result and increment the index of result for next edge
				if ( (x != y)  &&  e < (V - 1) )
				{
					Edge::vertices[minEdge.src].mst_degree++;
					Edge::vertices[minEdge.dest].mst_degree++;
					MST->edge[e++] = minEdge;
					Union(subsets, minEdge.x, minEdge.y);
				}
			}
		}

		v.push_back(next_edge);
		prevEdgeWeight = next_edge.weight;
	}

	free(subsets);

	return MST;
}


Graph* readDIMACS(string filename)
{
	string key1,key2;
	int vertices,edges;

	ifstream in(filename);
	in >> key1 >> key2;
	in >> vertices >> edges;
	Graph* graph = createGraph(vertices,edges);
	char e;
	int from, to, dist_val;
	int iter = 0;

	while(in >> e)
	{
		in >> from >> to >> dist_val;
		graph->edge[iter].src = from;
		graph->edge[iter].dest = to;
		graph->edge[iter].weight = dist_val;
		iter++;
	}

	return graph;
}

void dfs(vector<vector<int> > &adjList, bool* inBlock, int start)
{
	inBlock[start] = true;
	for (int i = 0; i < adjList.at(start).size(); ++i)
	{
		if (!inBlock[adjList.at(start).at(i)])
		{
			dfs(adjList, inBlock, adjList.at(start).at(i));
		}
	}
}

// Driver program to test above functions
int main( int argc, char *argv[] )
{
	if (argc != 3)
	{
		cout << "IMPROPER USAGE:  MUST HAVE [executable name] [filename of graph in DIMACS format] [desired number of blocks]" << endl;
		return 0;
	}

	// file containing graph in DIMACS format
	string filename;
	filename = argv[1];

	// figure out how many blocks we want
	int numPartitions;
	numPartitions = atoi(argv[2]);

	// build the graph object
	Graph* graph = readDIMACS( filename );

	// start the clock
	clock_t  start;
	start = clock();

	// find the MST of the graph
	Graph *MST = KruskalMST(graph, numPartitions);

	// convert to adjacency list
	vector<vector<int> > adjList, distances; // should be uninitialized at this point
	gimmeAdjListAndWeights(MST, adjList, distances);

	// find the rest of the distances
	vector<int> discovered;
	discovered.push_back(0);
	vector<int> drapers = adjList.at(0);
	findDist(adjList, distances, 0, drapers, discovered);

	// find the longest path
	vector<int> path = findPath(adjList, distances);

	int MSTweight = 0;
	for (int i = 0; i < MST->V - 1; ++i)
	{
		MSTweight += MST->edge[i].weight;
	}

	// run the dp algorithm
	long long ** DPcuts = partitionDPcuts(adjList, distances, path, MSTweight, numPartitions); // run the DP algorithm

	// recover the pieces
	long long  cutAfterIndex = path.size() - 1; // set cutAfterIndex to hold the info of the last edge in the path previously calculated
	vector <long long> bestWeights; // vector for storing the weights of the blocks of the graph
	vector <pair<int, int> > edgesOfSubPaths;
	pair<int, int> limitsOfCurrentBlock;
	vector<int>::iterator cutEdge;
	int left, right; // not necessary, but helps for readability
	while ( DPcuts[cutAfterIndex][1] != cutAfterIndex) // recover the edges to be cut
	{
		// RECALL:
		// theFavor[x][0] = running best penalty at edge x
		// theFavor[x][1] = the previous node after which a cut should be made which results in the penalty at theFavor[x][0] (best penalty achievable at this node)
		// theFavor[x][2] = total weight up to node x including branching off of node x but not including the weight between x and x + 1 (for calculating subtree weights)
		// theFavor[x][3] = weight of sub-portion of MST that occurs between nodes x and theFavor[x][1] (previous node that should be cut)
		// theFavor[x][4] = MSTweight - any edges (including edge x to x + 1) that have been cut to minimize penalty up to node x (for recalculating "ideal" block weight)

		bestWeights.push_back(DPcuts[cutAfterIndex][3]); // add the weight of the current block to the bestWeights vector

		limitsOfCurrentBlock.second = path.at(cutAfterIndex); // rightmost node in the subpath that constitutes this block
		limitsOfCurrentBlock.first = path.at(DPcuts[cutAfterIndex][1] + 1); // leftmost node in the subpath that constitutes this block
		// + 1 since DPcuts[cutAfterIndex][1] represents the previous cut-after-index which means the edge of this subpath is one index after that
		left = limitsOfCurrentBlock.first; // path.at(limitsOfCurrentBlock.first);
		right = limitsOfCurrentBlock.second; // path.at(limitsOfCurrentBlock.second);

		edgesOfSubPaths.push_back(limitsOfCurrentBlock); // store the pair that constitutes this block

		cutAfterIndex = DPcuts[cutAfterIndex][1]; // update cutAfterIndex to retrieve the next block at the next iteration of this loop

		// REMOVE THE EDGE FROM adjList
		// remove the right from the left list
		cutEdge = adjList.at(left).begin();
		for (int i = 0; i < adjList.at(left).size(); ++i, ++cutEdge)
		{
			if (*cutEdge == path.at(cutAfterIndex))
			{
				adjList.at(left).erase(cutEdge);
				break;
			}
		}
		// remove the left from the right list
		cutEdge = adjList.at( path.at(cutAfterIndex) ).begin();
		for (int i = 0; i < adjList.at( path.at(cutAfterIndex) ).size(); ++i, ++cutEdge)
		{
			if (*cutEdge == left)
			{
				adjList.at( path.at(cutAfterIndex) ).erase(cutEdge);
				break;
			}
		}
	}

	// push back the last block weight
	bestWeights.push_back(DPcuts[cutAfterIndex][3]);

	// build last set of limits on partial path...
	limitsOfCurrentBlock.first = path.at(0);
	limitsOfCurrentBlock.second = path.at(cutAfterIndex);

	// ... and push it onto the vector of these pairs
	edgesOfSubPaths.push_back(limitsOfCurrentBlock);
	left = limitsOfCurrentBlock.first;
	right = limitsOfCurrentBlock.second;

	// REMOVE THE LAST EDGE FROM adjList
	// remove the right from the left list
	if (right != path.back()) // account for the case when no cuts were made and we seek to return the whole MST... pointless, but necessary
	{
		cutEdge = adjList.at(right).begin();
		for (int i = 0; i < adjList.at(right).size(); ++i, ++cutEdge)
		{
			if (*cutEdge == path.at(cutAfterIndex + 1))
			{
				cout << "SHOULDN'T HAVE GOTTEN HERE" << endl;
				adjList.at(right).erase(cutEdge);
				break;
			}
		}
		// remove the left from the right list
		cutEdge = adjList.at( path.at(cutAfterIndex + 1) ).begin();
		for (int i = 0; i < adjList.at( path.at(cutAfterIndex + 1) ).size(); ++i, ++cutEdge)
		{
			if (*cutEdge == right)
			{
				cout << "SHOULDN'T HAVE GOTTEN HERE EITHER" << endl;
				adjList.at( path.at(cutAfterIndex + 1) ).erase(cutEdge);
				break;
			}
		}
	}

	// delete the DPcuts array
	for (int i = 0; i < path.size(); ++i)
	{
		delete[] DPcuts[i];
	}
	delete[] DPcuts;

	path.resize(0);


	// CUT MORE IF NECESSARY:
	int heaviestBlockIndex; // int for storing which block is the heaviest
	int infiniteLoopCatcher;
	vector<long long>::iterator weightIt, heaviestIt;
	vector<pair<int, int> >::iterator pairIt, heaviestPairIt;
	vector<int> partialPath(0); // the part of the path that constitutes the heaviest block (along which we want to make another cut)
	bool wasInWhile = false;

	while (bestWeights.size() < numPartitions) // if we don't have the required number of blocks...
	{
		wasInWhile = true;
		infiniteLoopCatcher = bestWeights.size();

		// ... find the heaviest block...
		heaviestBlockIndex = 0;
		partialPath.resize(0);
		weightIt = bestWeights.begin(); // reset weigtIt to the beginning of the bestWeights vector
		heaviestIt = bestWeights.begin(); // reset heaviestIt to the beginning of the bestWeights vector
		pairIt = edgesOfSubPaths.begin(); // reset pairIt to the beginning of the edgesOfSubPaths vector
		heaviestPairIt = edgesOfSubPaths.begin(); // reset heaviestPairIt to the beginning of the edgesOfSubPaths vector
		for (int i = 0; i < bestWeights.size(); ++i, ++weightIt, ++pairIt) // iterate over the blocks
		{
			if (bestWeights.at(i) > bestWeights.at(heaviestBlockIndex)) // check the weight of this block relative to running heaviest
			{
				heaviestBlockIndex = i; // this is now the running heaviest block
				heaviestIt = weightIt;
				heaviestPairIt = pairIt;
			}
		}

		// run dfs to find which nodes are in this block
		bool inBlock[adjList.size()]; // bool array for storing info about whether a node is in heaviest block
		for (int i = 0; i < adjList.size(); ++i)
		{
			inBlock[i] = false;
		}
		// now run dfs
		dfs(adjList, inBlock, (*heaviestPairIt).first);

		// build the part of the subpath that we want to cut again
		partialPath.resize(0);
		partialPath = miniFindPath(adjList, distances, inBlock);

		// ... run the DP algorithm on this block...
		long long ** cutPartialPath = partitionDPcuts(adjList, distances, partialPath, bestWeights.at(heaviestBlockIndex), 2);

		// erase the heaviest block in the current partition and the subpath in the edgesOfSubPaths vectors
		bestWeights.erase(heaviestIt);
		edgesOfSubPaths.erase(heaviestPairIt);

		// ... and update bestWeights and DPcuts with this info
		cutAfterIndex = partialPath.size() - 1;
		while ( cutPartialPath[cutAfterIndex][1] != cutAfterIndex )
		{
			// RECALL:
			// theFavor[x][0] = running best penalty at edge x
			// theFavor[x][1] = the previous node after which a cut should be made which results in the penalty at index 0 (best penalty achievable at this node)
			// theFavor[x][2] = total weight up to node x including branching off of node x but not including the weight between x and x + 1 (for calculating subtree weights)
			// theFavor[x][3] = weight of sub-portion of MST that occurs between nodes x and theFavor[x][1] (previous node that should be cut)
			// theFavor[x][4] = MSTweight - any edges (including edge x to x + 1) that have been cut to minimize penalty up to node x (for recalculating "ideal" block weight)

			bestWeights.push_back(cutPartialPath[cutAfterIndex][3]); // add the weight of the current block to the bestWeights vector

			limitsOfCurrentBlock.second = partialPath.at(cutAfterIndex); // rightmost node in the subpath that constitutes this block
			limitsOfCurrentBlock.first = partialPath.at(cutPartialPath[cutAfterIndex][1] + 1); // leftmost node in the subpath that constitutes this block
			// + 1 since DPcuts[cutAfterIndex][1] represents the previous cut-after-index which means the edge of this subpath is one index after that
			left = limitsOfCurrentBlock.first;
			right = limitsOfCurrentBlock.second;

			edgesOfSubPaths.push_back(limitsOfCurrentBlock); // store the pair that constitutes this block

			cutAfterIndex = cutPartialPath[cutAfterIndex][1]; // update cutAfterIndex to recrieve the next block at the next iteration of this loop

			// REMOVE THE EDGE FROM adjList
			// remove the right from the left list
			cutEdge = adjList.at(left).begin();
			for (int i = 0; i < adjList.at(left).size(); ++i, ++cutEdge)
			{
				if (*cutEdge == partialPath.at(cutAfterIndex))
				{
					adjList.at(left).erase(cutEdge);
					break;
				}
			}

			// remove the left from the right list
			cutEdge = adjList.at( partialPath.at(cutAfterIndex) ).begin();
			for (int i = 0; i < adjList.at( partialPath.at(cutAfterIndex) ).size(); ++i, ++cutEdge)
			{
				if (*cutEdge == left)
				{
					adjList.at( partialPath.at(cutAfterIndex) ).erase(cutEdge);
					break;
				}
			}
		}

		// push back the last block weight
		bestWeights.push_back(cutPartialPath[cutAfterIndex][3]);

		// build last set of limits on partial path...
		limitsOfCurrentBlock.first = partialPath.at(0);
		limitsOfCurrentBlock.second = partialPath.at(cutAfterIndex);

		// ... and push it onto the vector of these pairs
		edgesOfSubPaths.push_back(limitsOfCurrentBlock);
		left = limitsOfCurrentBlock.first;
		right = limitsOfCurrentBlock.second;


		// REMOVE THE LAST EDGE FROM adjList

		// remove the right from the left list
		if (right != partialPath.back()) // account for the case when no cuts were made and we seek to return the whole MST
		{
			cutEdge = adjList.at(right).begin();
			for (int i = 0; i < adjList.at(right).size(); ++i, ++cutEdge)
			{
				if (*cutEdge == partialPath.at(cutAfterIndex + 1))
				{
					adjList.at(right).erase(cutEdge);
					break;
				}
			}

			// remove the left from the right list
			cutEdge = adjList.at( partialPath.at(cutAfterIndex + 1) ).begin();
			for (int i = 0; i < adjList.at( partialPath.at(cutAfterIndex + 1) ).size(); ++i, ++cutEdge)
			{
				if (*cutEdge == right)
				{
					adjList.at(cutAfterIndex + 1).erase(cutEdge);
					break;
				}
			}
		}

		for (int i = 0; i < partialPath.size(); ++i)
		{
			delete[] cutPartialPath[i];
		}
		delete[] cutPartialPath;


		if (infiniteLoopCatcher == bestWeights.size())
		{
			cout << "POSSIBLE INFINITE LOOP DETECTED; BAILING OUT" << endl;
			exit(0);
		}
	}




	// PRINT OUT INFO

	// print the number of vertices and the number of edges
	cout << MST->V << "," << MST->E << ",";

	// print out block weights
	for (int i = 0; i < bestWeights.size() - 1; ++i)
	{
		cout << bestWeights.at(i) << "-";
	}
	cout << bestWeights.at(bestWeights.size() - 1); // print out the last block weight

	// print out the ratio of heaviest to lightest block weight
	cout <<","<<(double)*max_element(bestWeights.begin(),bestWeights.end()) / *min_element(bestWeights.begin(),bestWeights.end());


	// print out timing info
	cout<<","<<clock() - start / (double)(CLOCKS_PER_SEC)<<endl;

	free(graph);
	free(MST);

	return 0;
}