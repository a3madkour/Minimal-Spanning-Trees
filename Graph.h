#ifndef __GRAPH__
#define __GRAPH__
#include <string>
#include <fstream>
#include <vector>

using namespace std;

struct Vertex
{
    int id;
    int mst_degree = 0;
};

struct subset
{
    int parent;
    int rank;
};
struct Edge
{

    static vector<Vertex> vertices;

    int src, dest, weight, x,y;
    bool operator==(const Edge& rhs)
{
    return (src==rhs.src && dest==rhs.dest && weight==rhs.weight);
}
   bool operator>=(const Edge& rhs)
{
    return (weight>=rhs.weight);
}
   bool operator>(const Edge& rhs)
{
    return (weight>rhs.weight);
}
   bool operator<(const Edge& rhs) const
{
    return (weight<rhs.weight);
}
 Edge(){
        
    }
};

class Graph
{
private:
	// V-> Number of vertices, E-> Number of edges
    int V, E;
    // graph is represented as an array of edges. Since the graph is
    // undirected, the edge from src to dest is also edge from dest
    // to src. Both are counted as 1 edge here.
    vector <Edge> MST;
    Graph * parent;
    int MSTweight;
    int find(struct subset subsets[], int i);
    void Union(struct subset subsets[], int x, int y);
    static int myComp(const void* a, const void* b);
    int minDistance(int dist[], bool sptSet[], int vertices);
    void findPath (int parent[],int j,int src,vector <int> &path);
protected:
  int ** adjMat;
public:
    struct Edge* edge;
    void dijkstra(int src, vector<Edge> &newEdges, bool useParent=false);
    Graph(int v, int e, Graph * Parent);
    Graph(string filename);
    Graph();
    int * indices;
    void assignEdges(int e);
    void setIndex (int original, int altered);
    void toAdjMat ();
    void setParent(Graph * Parent){parent = Parent;}
    void printAdjMat(); 
    void printEdges();
    int getV(){return V;}
    int getE(){return E;}
    Edge * getEdge(){return edge;}
    vector <Edge> getMST(){return MST;}
    ~Graph();
    int modKruskalMST();
    void readSegments(Graph *negPartition, Graph *posPartition, double * eigenVector);
};

#endif // __GRAPH__ 
