#include "Graph.h"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <limits.h>

using namespace std;
ostream& operator<<(ostream& os, Edge& e )
{
    os << e.src << "---" << e.dest << " == " << e.weight;
    return os;
}

vector<Vertex> Edge :: vertices;
Graph:: Graph (int v, int e, Graph* Parent=0){

    V = v;
    E = e;
    adjMat = 0;
    parent = Parent;
    for (int i = 0; i < V; ++i)
    {
        Vertex v;
        v.id = i;
        Edge::vertices.push_back(v);
    }
    indices = new int [V];

}

Graph :: Graph(string filename){
    string key1,key2;
    int vertices,edges;

    ifstream in(filename);
    in >> key1 >> key2;
    in >> vertices >> edges;
    parent = 0;
    V = vertices;
    E = edges;
    adjMat = 0;
    indices = new int [V];
    for (int i = 0; i < V; ++i)
    {
        indices[i] = i;
        Vertex v;
        v.id = i;
        Edge::vertices.push_back(v);
    }
    edge = new Edge [E];

    char e;
    int from, to, dist_val;
    int iter = 0;
    while(in >> e)
    {
        in >> from >> to >> dist_val;
        edge[iter].src = from;
        edge[iter].dest = to;
        edge[iter].weight = dist_val;
        iter++;

    }
     toAdjMat();
}

Graph :: Graph(){
    adjMat = 0;
    parent = 0;
    edge = 0;
    indices = 0;
}
void Graph :: toAdjMat(){
    adjMat = new int*[V];
    for (int i = 0; i < V; ++i)
    {
        adjMat[i] = new int[V];
    }
   for (int i = 0; i < V; ++i)
   {
       for (int j = 0; j < V; ++j)
       {
           adjMat[i][j] = 0;
       }
   }
    for (int i = 0; i < E; ++i)
    {
        adjMat[edge[i].src][edge[i].dest] = edge[i].weight;
        adjMat[edge[i].dest][edge[i].src] = edge[i].weight;
        adjMat[edge[i].src][edge[i].src] = 0;
    }

}
void Graph :: printAdjMat(){
    if(adjMat){
        for (int i = 0; i < V; ++i)
        {
            for (int j = 0; j < V; ++j)
            {
                cout<<adjMat[i][j]<<" ";
            }
            cout<<endl;
        }
    }else{
        cout<<"No adjMat"<<endl;
    }


}
void Graph :: readSegments(Graph *negPartition, Graph *posPartition, double * eigenVector)
{


    vector<bool> inNegativePartition;
    negPartition->setParent(this);
    posPartition->setParent(this);
    int numberVertsInNegPartition = 0;
    for (int i = 0; i < V; ++i)
    {
        if (eigenVector[i] < 0)
        {
            inNegativePartition.push_back(true);

            ++numberVertsInNegPartition;
        }
        else
        {
            inNegativePartition.push_back(false);
        }
    }


    // figure out how many edges will go in each partition
    int numEdgesNegPartition = 0, numEdgesPosPartition = 0;
    for (int i = 0; i < E; ++i)
    {
        if ( inNegativePartition.at(edge[i].src) && inNegativePartition.at(edge[i].dest) )
        {
            numEdgesNegPartition++;
        }
        else if ( !inNegativePartition.at(edge[i].src) && !inNegativePartition.at(edge[i].dest) )
        {
            numEdgesPosPartition++;
        }
    }

    // populate the edge lists of each partition
    negPartition->V = numberVertsInNegPartition, posPartition->V = V - numberVertsInNegPartition;
    negPartition->E = numEdgesNegPartition, posPartition->E = numEdgesPosPartition;
    negPartition->edge = new Edge[numEdgesNegPartition];
    posPartition->edge = new Edge[numEdgesPosPartition];
    int negPartCount = 0, posPartCount = 0;
    for (int i = 0; i < E; ++i)
    {
        if ( inNegativePartition.at(edge[i].src) && inNegativePartition.at(edge[i].dest) )
        {
            negPartition->edge[negPartCount] = edge[i];
            ++negPartCount;
        }
        else if ( !inNegativePartition.at(edge[i].src) && !inNegativePartition.at(edge[i].dest) )
        {
            posPartition->edge[posPartCount] = edge[i];
            ++posPartCount;
        }
    }

    //store the indicies
    int newNegIndex= 0;
    int newPosIndex = 0;
    for (int i = 0; i < V; ++i)
    {
        if (eigenVector[i] < 0)
        {

            negPartition->setIndex(i,newNegIndex);
            newNegIndex++;

        }
        else
        {
            posPartition->setIndex(i,newPosIndex);
            newPosIndex++;
        }
    }
    // swap out the src and dest to correspond to 0-based sequential indices instead of random all over the place
    int negIndex = 0, posIndex = 0;
    vector<int> newPartNums(inNegativePartition.size());
    for (int i = 0; i < inNegativePartition.size(); ++i)
    {
        if (inNegativePartition.at(i))
        {
            newPartNums.at(i) = negIndex++;
        }
        else
        {
            newPartNums.at(i) = posIndex++;
        }
    }
    for (int i = 0; i < negPartition->E; ++i)
    {

        negPartition->edge[i].src = newPartNums.at(negPartition->edge[i].src);
        negPartition->edge[i].dest = newPartNums.at(negPartition->edge[i].dest);
    }
    for (int i = 0; i < posPartition->E; ++i)
    {

        posPartition->edge[i].src = newPartNums.at(posPartition->edge[i].src);
        posPartition->edge[i].dest = newPartNums.at(posPartition->edge[i].dest);
    }

}
void Graph :: printEdges(){
    cout<<"Edges: "<<endl;
    for (int i = 0; i < E; ++i)
    {
        cout<<edge[i]<<endl;
    }

}
Graph :: ~Graph(){
    if(edge){
        delete [] edge;
    }
    if(adjMat){
        for (int i = 0; i < V; ++i)
        {
            delete [] adjMat[i];
        }
        delete [] adjMat;
    }
    if(indices){
        delete [] indices;
    }
}

void Graph :: setIndex (int original, int altered){
    if(indices==0){
        indices = new int [V];
    }
    indices[altered] = original;

}

// A utility function to find set of an element i
// (uses path compression technique)
int Graph :: find(struct subset subsets[], int i)
{
    // find root and make root as parent of i (path compression)
    if (subsets[i].parent != i)
        subsets[i].parent = find(subsets, subsets[i].parent);

    return subsets[i].parent;
}


void Graph :: assignEdges(int e){

	if(edge){
		delete edge;
	}

	edge = new Edge [e];

}

// A function that does union of two sets of x and y
// (uses union by rank)
void Graph :: Union(struct subset subsets[], int x, int y)
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
int Graph :: myComp(const void* a, const void* b)
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


// The main function to construct MST using Kruskal's algorithm
int Graph :: modKruskalMST()
{

    int e = 0; // An index variable, used for result[]
    int i = 0; // An index variable, used for sorted edges
    MSTweight = 0;
    // Step 1: Sort all the edges in non-decreasing order of their weight
    // If we are not allowed to change the given graph, we can create a copy of
    // array of edges
    qsort(edge, E, sizeof(edge[0]), myComp);

    // Allocate memory for creating V ssubsets
    struct subset *subsets =
        (struct subset*) malloc( V * sizeof(struct subset) );

    // Create V subsets with single elements
    for (int v = 0; v < V; ++v)
    {
        subsets[v].parent = v;
        subsets[v].rank = 0;
    }
    std::vector<Edge> edgeVector;
    int prevEdgeWeight = -999999;
    // Number of edges to be taken is equal to V-1
    while (e < V - 1)
    {
        // Step 2: Pick the smallest edge. And increment the index
        // for next iteration
        if(i > E){
            for (int j = 0; j < V; ++j)
            {
                if(subsets[j].parent == j && subsets[j].rank == 0){
                    vector <Edge> newEdges;
                    dijkstra(j,newEdges,true);
                    for (int i = 0; i < newEdges.size(); ++i)
                    {

                        MSTweight+= newEdges[i].weight;
                        MST.push_back(newEdges[i]);
                    }

                }
            }
            return MSTweight;
        }
        struct Edge next_edge = edge[i++];

        if(prevEdgeWeight!= next_edge.weight && prevEdgeWeight > 0){
            int minDegree;
            int currDegree;
            int index;
            Edge minEdge;
            minEdge.src = -1;
            while(edgeVector.size()!=0){
                minDegree = 31231231;
                for (int j = 0; j < edgeVector.size(); ++j)
                {

                    currDegree = Edge::vertices[edgeVector[j].src].mst_degree + Edge::vertices[edgeVector[j].dest].mst_degree;
                    if(currDegree < minDegree){
                        minDegree = currDegree;
                        minEdge = edgeVector[j];
                        index = j;
                        //cout<<"minDegree: "<<minDegree<<endl;
                        //cout<<"minEdge: "<<minEdge.src<< "--->" <<  minEdge.dest<< endl;
                    }

                }
                edgeVector.erase(edgeVector.begin()+index);
                int x = find(subsets, minEdge.src);
                int y = find(subsets, minEdge.dest);
                minEdge.x = x;
                minEdge.y = y;
                // If including this edge does't cause cycle, include it
                // in result and increment the index of result for next edge
                if (x != y  &&  e < V - 1)
                    {
                    Edge::vertices[minEdge.src].mst_degree++;
                    Edge::vertices[minEdge.dest].mst_degree++;
                    MST.push_back(minEdge);
                    e++;
                    MSTweight+=minEdge.weight;
                    Union(subsets, minEdge.x, minEdge.y);
                    }

            }
        }
        edgeVector.push_back(next_edge);
        prevEdgeWeight = next_edge.weight;
    }

    // print the contents of result[] to display the built MST
    // cout << "graph->V = " << graph->V << endl;
    // printf("Following are the edges in the constructed MST\n");
    // for (i = 0; i < graph->V - 1; ++i)
    //     printf("%d -- %d == %d\n", result[i].src, result[i].dest,
    //                                             result[i].weight);

    return MSTweight;
}

int Graph :: minDistance(int dist[], bool sptSet[], int vertices)
{
    // Initialize min value
    int min = INT_MAX, min_index;

    for (int v = 0; v < vertices; v++)
        if (sptSet[v] == false && dist[v] <= min)
            min = dist[v], min_index = v;

    return min_index;
}

void Graph :: findPath (int parent[],int j,int src,vector <int> &path){
    //find edge using adjmat
    // Base Case : If j is source
    if (parent[j]==-1)
        return;
    findPath(parent, parent[j],src,path);

    path.push_back(j);
}

void Graph ::dijkstra(int src, vector <Edge> &newEdges, bool useParent){

    Graph * graph;
    int indexSrc = src;
    if(useParent){
        if(parent == 0){
                //cout<<"Graph has no parent"<<endl;
                return ;
        }
        if(getE()==0 || getV()==0 || indices == 0){
                //cout<<"Something went wrong"<<endl;
                return;
        }
        graph = parent;
        while(graph->parent!=0){
            graph = graph->parent;
            for (int i = 0; i < V; ++i)
            {
                indices[i] = graph->indices[indices[i]];
            }
        }
        indexSrc = indices[src];
    }else{
        graph = this;
    }


    int dist[graph->getV()];  // The output array. dist[i] will hold
                  // the shortest distance from src to i

    // sptSet[i] will true if vertex i is included / in shortest
    // path tree or shortest distance from src to i is finalized
    bool sptSet[graph->getV()];

    // Parent array to store shortest path tree
    int parent[graph->getV()];

    // Initialize all distances as INFINITE and stpSet[] as false
    for (int i = 0; i < graph->getV(); i++)
    {
        parent[indexSrc] = -1;
        dist[i] = INT_MAX;
        sptSet[i] = false;
    }

    // Distance of source vertex from itself is always 0
    if(useParent){
         dist[indexSrc] = 0;
     }else{
         dist[src] = 0;
     }


    // Find shortest path for all vertices

    for (int count = 0; count < graph->getV()-1; count++)
    {
        // Pick the minimum distance vertex from the set of
        // vertices not yet processed. u is always equal to src
        // in first iteration.
        int u = minDistance(dist, sptSet, graph->getV());

        // Mark the picked vertex as processed
        sptSet[u] = true;

        // Update dist value of the adjacent vertices of the
        // picked vertex.
        for (int v = 0; v < graph->getV(); v++){

            // Update dist[v] only if is not in sptSet, there is
            // an edge from u to v, and total weight of path from
            // src to v through u is smaller than current value of
            // dist[v]

            if (!sptSet[v] && graph->adjMat[u][v] &&
                dist[u] + graph->adjMat[u][v] < dist[v])
            {

                parent[v]  = u;
                dist[v] = dist[u] + graph->adjMat[u][v];

            }
        }
    }
    int dest = src;
    if(useParent){
        int minDist = INT_MAX;
        for (int i = 0; i < graph->getV(); ++i)
        {

            for (int j = 0; j <V; ++j)
            {
                if(i == indices[j] && i!=src){
                    if(dist[i]< minDist){

                            minDist = dist[i];
                            dest = i;
                    }
                }
            }
        }

        vector <int> path;
        path.push_back(src);
        findPath(parent,dest,indexSrc,path);
        for (int i = 0; i < path.size()-1; ++i)
        {
            Edge e;
            e.src= path[i];
            e.dest = path[i+1];
            e.weight = graph->adjMat[e.src][e.dest];
            e.src += graph->getV();
            e.dest += graph->getV();
            newEdges.push_back(e);
        }
    }


    // print the constructed distance array
    //printf("Vertex\t  Distance\tPath");
    for (int i = 1; i < graph->getV(); i++)
    {
        vector <int> path;

     //   printf("\n%d -> %d \t\t %d\t\t%d ", src, i, dist[i], src);
        findPath(parent, i,src,path);

    }
    //cout<<endl;

    return;
}
