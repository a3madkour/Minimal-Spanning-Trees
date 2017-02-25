
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <ctime>
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
	Edge& operator= (const Edge& e){
		src = e.src;
		dest = e.dest;
		weight = e.weight;
		x = e.x;
		y = e.y;
	}

};

ostream& operator<<(ostream& os, Edge& e )
{
    os << e.src << "---" << e.dest << " == " << e.weight;
    return os;
}

vector<Vertex> Edge :: vertices;

struct Cuts
{
	vector <Edge> edges;
	vector <int>  partitions;
	bool operator<(const Edge& rhs)
	{
    	bool smaller;
    	return smaller;
	}
	Cuts replace(Edge e, int index){
		edges[index] = e;
		find_paritions();
	}
	void find_paritions(){
		cout<<"Finding partitions"<<endl;
	}
	Cuts& operator= (const Cuts& cuts){
		edges.resize(cuts.edges.size());
		partitions.resize(cuts.partitions.size());
		for (int i = 0; i < edges.size(); ++i)
		{
			edges[i] = cuts.edges[i];
		}
		for (int i = 0; i < partitions.size(); ++i)
		{
			partitions[i] = cuts.partitions[i];
		}
	}
	Cuts(int k){
		edges.resize(k-1);
		partitions.resize(k);
	}
	Cuts(vector <Edge> Edges, vector <int> Partitions){
		edges.resize(Edges.size());
		partitions.resize(Partitions.size());
		for (int i = 0; i < edges.size(); ++i)
		{
			edges[i] = Edges[i];
		}
		for (int i = 0; i < partitions.size(); ++i)
		{
			partitions[i] = Partitions[i];
		}
	}
	Cuts(const Cuts& cuts){
		edges.resize(cuts.edges.size());
		partitions.resize(cuts.partitions.size());
		for (int i = 0; i < edges.size(); ++i)
		{
			edges[i] = cuts.edges[i];
		}
		for (int i = 0; i < partitions.size(); ++i)
		{
			partitions[i] = cuts.partitions[i];
		}
	}


};

ostream& operator<<(ostream& os, Cuts& cuts )
{
    os << "Edges Cut: "<<endl;
    for (int i = 0; i < cuts.edges.size(); ++i)
    {
    	os<<cuts.edges[i]<<endl;
    }
    os << "Weights of the partitions: "<<endl;
    for (int i = 0; i < cuts.partitions.size()-1; ++i)
    {
    	os<<cuts.partitions[i]<<"-";
    }
    os<<cuts.partitions[cuts.partitions.size()-1]<<endl;
    return os;
}

// a structure to represent a connected, undirected and weighted graph
struct Graph
{
	// V-> Number of vertices, E-> Number of edges
	int V, E;

	// graph is represented as an array of edges. Since the graph is
	// undirected, the edge from src to dest is also edge from dest
	// to src. Both are counted as 1 edge here.
	struct Edge* edge;
};
// Creates a graph with V vertices and E edges
struct Graph* createGraph(int V, int E)
{
	struct Graph* graph = (struct Graph*) malloc( sizeof(struct Graph) );
	graph->V = V;
	graph->E = E;

	for (int i = 0; i < V; ++i)
	{
		Vertex v;
		v.id = i;
		Edge::vertices.push_back(v);
	}
	graph->edge = (struct Edge*) malloc( graph->E * sizeof( struct Edge ) );

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
int find(struct subset subsets[], int i)
{
	// find root and make root as parent of i (path compression)
	if (subsets[i].parent != i)
		subsets[i].parent = find(subsets, subsets[i].parent);

	return subsets[i].parent;
}


// A function that does union of two sets of x and y
// (uses union by rank)
void Union(struct subset subsets[], int x, int y)
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

void  KruskalMST(struct Graph* graph, vector <Edge> &result)
{
	int V = graph->V;
	int e = 0; // An index variable, used for result[]
	int i = 0; // An index variable, used for sorted edges

	// Step 1: Sort all the edges in non-decreasing order of their weight
	// If we are not allowed to change the given graph, we can create a copy of
	// array of edges
	qsort(graph->edge, graph->E, sizeof(graph->edge[0]), myComp);

	// Allocate memory for creating V ssubsets
	struct subset *subsets = (struct subset*) malloc( V * sizeof(struct subset) );

	// Create V subsets with single elements
	for (int v = 0; v < V; ++v)
	{
		subsets[v].parent = v;
		subsets[v].rank = 0;
	}
	std::vector<Edge> v;
	int prevEdgeWeight = -999999;
	// Number of edges to be taken is equal to V-1
	while (e < V - 1)
	{
		// Step 2: Pick the smallest edge. And increment the index
		// for next iteration
		struct Edge next_edge = graph->edge[i++];

		if(prevEdgeWeight!= next_edge.weight && prevEdgeWeight > 0)
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

					currDegree = Edge::vertices[v[j].src].mst_degree + Edge::vertices[v[j].dest].mst_degree;
					if(currDegree < minDegree)
					{
						minDegree = currDegree;
						minEdge = v[j];
						index = j;
						//cout<<"minDegree: "<<minDegree<<endl;
						//cout<<"minEdge: "<<minEdge.src<< "--->" <<  minEdge.dest<< endl;
					}

				}
				v.erase(v.begin()+index);
				//cout<<"v.size(): "<<v.size()<<endl;
				int x = find(subsets, minEdge.src);
				int y = find(subsets, minEdge.dest);
				minEdge.x = x;
				minEdge.y = y;
				// If including this edge does't cause cycle, include it
				// in result and increment the index of result for next edge
				//cout<<"mamamia 3"<<endl;
				if (x != y  &&  e < V - 1)
				{
					Edge::vertices[minEdge.src].mst_degree++;
					Edge::vertices[minEdge.dest].mst_degree++;
					//cout<<"mamamia 2"<<endl;
					result[e++] = minEdge;
					Union(subsets, minEdge.x, minEdge.y);
				}
			}
		}

		v.push_back(next_edge);
		prevEdgeWeight = next_edge.weight;
	}


	// print the contents of result[] to display the built MST
	// printf("Following are the edges in the constructed MST\n");
	// for (i = 0; i < e; ++i)
	// 	printf("%d -- %d == %d\n", result[i].src, result[i].dest, result[i].weight);

	return /*result*/;
}
struct Graph* readDIMACS(string filename)
{
	string key1,key2;
	int vertices,edges;

	ifstream in(filename);
	in >> key1 >> key2;
	in >> vertices >> edges;
	struct Graph* graph = createGraph(vertices,edges);
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
int main( int argc, char *argv[] )
{
	 int num_partitions;
     struct Graph* graph;
     clock_t  start;
     if(argc == 1){
           cout<<"No graph specified"<<endl;
           // return 0;
     }else if(argc == 2){
           graph = readDIMACS(argv[1]);
           num_partitions = 3;
     }else{
           stringstream convert(argv[2]);
           graph = readDIMACS(argv[1]);
           convert >> num_partitions;
     }



}