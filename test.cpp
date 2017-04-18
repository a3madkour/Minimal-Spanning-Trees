/* system example : DIR */
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* system, NULL, EXIT_FAILURE */
#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include "prettyPrint.h"

using namespace std;
int main ()
{
	int i;
	printf ("Checking if processor is available...");
	int max_num_partitions  = 8;
	if (system(NULL)) puts ("Ok");
	else exit (EXIT_FAILURE);
	vector<char> sample;
	sample.push_back('A');
	sample.push_back('B');
	sample.push_back('C');
	string filename,filename2;
	string output ("output.csv");
	string command;
	clock_t  start;
	system("echo  File Name,Algorithim,Number of Partitions,Vertices, Edges,Weights of Partitions,Weight of MST,Multiplicative Spread,Elapsed Time >> testData.csv");
	for (int i = 40; i <= 100; i+=10)
	{
		//do stuff "graph"+i+"Afull.txt"
		for (int j = 0; j < sample.size(); ++j)
		{
			for (int k = 2; k <= max_num_partitions; k++)
			{

			filename = "graph"  +to_string(i) + sample[j] + "half.txt";
			command = "echo -n "+filename+",Dynamic Programming Approach with All Trees," + to_string(k)+", >> "+output;
			system(command.c_str());
			command  = "./newdp graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			system(command.c_str());

			filename = "graph"  +to_string(i) + sample[j] + "full.txt";
			command = "echo -n "+filename+",Dynamic Programming Approach with All Trees," + to_string(k)+", >> "+output;
			system(command.c_str());
			command  = "./newdp graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			system(command.c_str());

			filename = "graph"  +to_string(i) + sample[j] + "half.txt";
			command = "echo -n "+filename+",Segmantation Based Clustring with Modified Kruskal," + to_string(k)+", >> "+output;
			system(command.c_str());
			command  = "./ncut graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			system(command.c_str());

			cout<<"We are at: "<<filename<<endl;

			filename = "graph"  +to_string(i) + sample[j] + "full.txt";
			command = "echo -n "+filename+",Segmantation Based Clustring with Modified Kruskal," + to_string(k)+", >> "+output;
			system(command.c_str());
			command  = "./ncut graphs/"+filename + " " +to_string(k)+ " >> "+ output;
			system(command.c_str());
			}
		}
	}
	return 0;
}
