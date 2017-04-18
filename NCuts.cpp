
#include <lapacke.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <ctime>
#include <numeric>
#include "Graph.h"


using namespace std;


void printMatrix(vector <vector <int> > &matrix)
{
    for (int i = 0; i < matrix.size(); ++i)
    {
        for (int j = 0; j < matrix[i].size(); ++j)
        {
            cout<<matrix[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<< endl;
    return;
}

void floydWarshall(Graph * graph, vector < vector < int > > & dist){
    for (int i = 0; i < graph->getV(); ++i)
    {
        vector < int > v;
        for (int j = 0; j < graph->getV(); ++j)
        {
            v.push_back(99999999);

        }
        dist.push_back(v);
        dist[i][i] = 0;
    }

    Edge * graphEdge = graph->getEdge();
    for (int i = 0; i < graph->getE(); ++i)
    {
        dist[graphEdge[i].src][graphEdge[i].dest] = graphEdge[i].weight;
        dist[graphEdge[i].dest][graphEdge[i].src] = graphEdge[i].weight;
    }

    for (int k = 0; k < graph->getV(); ++k)
    {
        for (int i = 0; i < graph->getV(); ++i)
        {
            for (int j = 0; j < graph->getV(); ++j)
            {
                if(dist[i][j] > dist[i][k] + dist[k][j])
                {


                    dist[i][j] = dist[i][k] + dist[k][j];

                }
            }
        }
    }
}

void calculateW(vector < vector < int> > & W, vector < vector < int> > & dist,Graph * graph){
    int max_weight = -1;
    for (int i = 0; i < dist.size(); ++i)
    {
        int temp = (*max_element(dist[i].begin(),dist[i].end()));
        if(temp > max_weight){
                max_weight = temp;
        }
    }
    max_weight++;
    for (int i = 0; i < graph->getV(); ++i)
    {
        for (int j = 0; j < graph->getV(); ++j)
        {
            if(i==j){
                W[i][j] = 0;
            }else{
                W[i][j] = max_weight - dist[i][j];
            }

        }

    }
}
void calculateD(vector < vector < int> > & D, vector < vector  < int > > & W,vector < vector < int> > & dist, Graph * graph){

    for (int i = 0; i < graph->getV(); ++i)
    {


            D[i][i] = accumulate(W[i].begin(),W[i].end(),0);

    }
}

bool isPowerOfTwo(int n)
{
  if (n == 0)
    return 0;
  while (n != 1)
  {
    if (n%2 != 0)
      return 0;
    n = n/2;
  }
  return 1;
}

void NCut (Graph * graph, int num_partitions, vector <int>& weights){
     clock_t start;
     start = clock();
     vector < vector < int > > dist;
     floydWarshall(graph,dist);
     clock_t endFlyod = clock()-start;
     cout<<"End of floyd-warshall: " << (clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" <<endl;   
     vector < vector < int> > W ;
     vector < vector < int> > D ;
     vector < vector < int> > Ma ;

     for (int i = 0; i < graph->getV(); ++i)
     {
        vector <int> v;
        for (int j = 0; j < graph->getV(); ++j)
             {
                v.push_back(0);
             }

              W.push_back(v);
              D.push_back(v);
              Ma.push_back(v);
     }
    calculateW(W,dist,graph);
    calculateD(D,W,dist,graph);
    // note, to understand this part take a look in the MAN pages, at section of parameters.

    int     N = graph->getV();
    int     ITYPE = 1;
    char    JOBZ = 'V';
    char    RANGE = 'I';
    char    UPLO  = 'U';
    int     INFO=3;
    double  VL = 2;
    double  VU = 2;
    int     IL = 2;
    int     IU = 2;
    double  ABSTOL =  0.01;
    int     M = 1;
    double * Wi = new double [N];
    double  * Z = new double [N];
    int     LDZ = N;
    int     * IFAIL = new int [N];




    double * AP  = new double [N*(N+1)/2];
    double * BP  = new double [N*(N+1)/2];
// end of declarations
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            Ma[i][j] = D[i][j]-W[i][j];
        }
    }

    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i <= j; ++i)
        {

            AP[i + j*(j+1)/2] = Ma[i][j];
            BP[i + j*(j+1)/2] = D[i][j];
        }

    }

    INFO = LAPACKE_dspgvx(LAPACK_COL_MAJOR,ITYPE,JOBZ,RANGE,UPLO,N,AP,BP,VL,VU,IL,IU,ABSTOL,&M,Wi,Z,LDZ,IFAIL);
    if(INFO)
    {
    }else{
    }

    Graph negPartition, posPartition;
    graph->readSegments(&negPartition, &posPartition, Z);


    Edge * posPartitionEdge = posPartition.getEdge();

    Edge * negPartitionEdge = negPartition.getEdge();
    int posWeight = posPartition.modKruskalMST();
    int negWeight = negPartition.modKruskalMST();

    cout<<"End of ncut: " << (clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" <<endl;   

    if(num_partitions == 2){
        weights.push_back(posWeight);
        weights.push_back(negWeight);   
        return;
    }else{
        if(posPartition.getE()<2){
            weights.push_back(posWeight);
            cout<<"Can't subdivide further"<<endl;
        }else{
            NCut(&posPartition,num_partitions/2,weights);
        }
        
        if (negPartition.getE()<2)
        {
            weights.push_back(negWeight);
            cout<<"Can't subdivide further"<<endl;
        }else{
            NCut(&negPartition,num_partitions/2,weights);
        }
        return;
    }



}
int main(int argc, char * argv[])
{
    Graph * graph;
     int num_partitions = 2;
     if(argc == 1){
           cout<<"No graph specified"<<endl;
            graph = new Graph("graphs/graph10Afull.txt");
            num_partitions = 16;
     }else if(argc == 2){
                graph = new Graph(argv[1]);


     }else if(argc == 3){
            graph = new Graph(argv[1]);
            if(isPowerOfTwo(atoi(argv[2]))){
                    num_partitions = atoi(argv[2]);
             }


     }

     int weight = graph->modKruskalMST();
     vector <int> weights;
     NCut(graph,num_partitions,weights);
     for (int i = 0; i < weights.size()-1; ++i)
     {
         cout<<weights[i]<<"-";
     }
     cout<<weights[weights.size()-1]<<",";
     cout<< (double) *max_element(weights.begin(),weights.end()) /weight<<",";
     cout<<weight<<",";
     cout<< weight /num_partitions <<endl;

    return 0;
}
