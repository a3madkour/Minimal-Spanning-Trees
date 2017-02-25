#include <iostream>	
#include "Graph.h"

using namespace std;

int main(){
	
	Graph graph ("graph351.txt");
	cout<<"MSTweight: "<<graph.modKruskalMST()<<endl;
	std::vector<Edge> v = graph.getMST();
	// for (int i = 0; i < v.size(); ++i)
	// {
	// 	cout<<v[i].src<<" "<<v[i].dest<<endl;
	// }
	double Z[] = {-1,1,1,1,1,-1,-1,-1,-1};
	Graph negPartition, posPartition;
    graph.readSegments(&negPartition, &posPartition, Z);
    cout<<"neg"<<endl;
    cout<<negPartition.getE()<<endl;
	negPartition.toAdjMat();
	vector <Edge> test;
	cout<<"V: "<<negPartition.getV()<<endl;
	for (int i = 0; i < negPartition.getV(); ++i)
	{
		cout<<negPartition.indices[i]<<" ";
	}
	negPartition.printEdges();
	cout<<endl;
	negPartition.dijkstra(0,test,true);
	// graph.printAdjMat();
	cout<<"Test: "<<endl;
	for (int i = 0; i < test.size(); ++i)
	{
		cout<<test[i].src<<"   "<<test[i].dest<<" === "<<test[i].weight<<endl;
	}
	return 0;
}