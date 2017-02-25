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
	string output ("kaptin11.csv");
	string command;
	clock_t  start;
	system("echo  File Name,Algorithim,Number of Partitions,Vertices, Edges,Weights of Partitions,Weight of MST,Multiplicative Spread,Elapsed Time >> testData.csv");
	// your test

	//   for (int i = 2; i <= 7; i++)
	// {
	// 		//do stuff "graph"+i+"Afull.txt"
	// 		for (int j = 0; j < sample.size(); ++j)
	// 		{
	// 			for (int k = 2; k <= max_num_partitions; ++k)
	// 			{

	// 			filename = "graph"  +to_string(i) + sample[j] + "half.txt";
	// 			command = "echo -n "+filename+",Brute Force," + to_string(k)+", >> "+output;
	// 			system(command.c_str());
	// 			command  = "./bruteForce2 graphs/"+filename + " " +to_string(k)+ " >> "+ output;
	// 			system(command.c_str());

	// 			filename = "graph"  +to_string(i) + sample[j] + "full.txt";
	// 			command = "echo -n "+filename+",Brute Force," + to_string(k)+", >> "+output;
	// 			system(command.c_str());
	// 			command  = "./bruteForce2 graphs/"+filename + " " +to_string(k)+ " >> "+ output;

	// 			system(command.c_str());


	// 			}
	// 		}
	// }
	for (int i = 40; i <= 100; i+=10)
	{
		//do stuff "graph"+i+"Afull.txt"
		for (int j = 0; j < sample.size(); ++j)
		{
			for (int k = 2; k <= max_num_partitions; k++)
			{

			// filename = "graph"  +to_string(i) + sample[j] + "half.txt";
			// command = "echo -n "+filename+",Consecutive Cuts with Modified Kruskal," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./ConsecutiveCutsModKrus graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			// system(command.c_str());

			// cout<<"We are at: "<<filename<<endl;

			// filename = "graph"  +to_string(i) + sample[j] + "full.txt";
			// command = "echo -n "+filename+",Consecutive Cuts with Modified Kruskal," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./ConsecutiveCutsModKrus graphs/"+filename + " " +to_string(k)+ " >> "+ output;
			// system(command.c_str());

			// filename = "graph"  +to_string(i) + sample[j] + "half.txt";
			// command = "echo -n "+filename+",Consecutive Cuts with All trees," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./ConsecutiveCutsAllTrees graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			// system(command.c_str());

			// filename = "graph"  +to_string(i) + sample[j] + "full.txt";
			// command = "echo -n "+filename+",Consecutive Cuts with All trees," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./ConsecutiveCutsAllTrees graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			// system(command.c_str());



			// filename = "graph"  +to_string(i) + sample[j] + "half.txt";
			// command = "echo -n "+filename+",Recursive Cuts with Modified Kruskal," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./RecursiveCutsModKrus graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			// system(command.c_str());

			// filename = "graph"  +to_string(i) + sample[j] + "full.txt";
			// command = "echo -n "+filename+",Recursive Cuts with Modified Kruskal," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./RecursiveCutsModKrus graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			// system(command.c_str());

			//  filename = "graph"  +to_string(i) + sample[j] + "half.txt";
			// command = "echo -n "+filename+",Recursive Cuts with All trees," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./RecursiveCutsAllTrees graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			// system(command.c_str());

			// filename = "graph"  +to_string(i) + sample[j] + "full.txt";
			// command = "echo -n "+filename+",Recursive Cuts with All trees," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./RecursiveCutsAllTrees graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			// system(command.c_str());



			// filename = "graph"  +to_string(i) + sample[j] + "half.txt";
			// command = "echo -n "+filename+",Dynamic Programming Approach with Modified Kruskal," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./DynamicProgrammingCutsModKrus graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			// system(command.c_str());

			// filename = "graph"  +to_string(i) + sample[j] + "full.txt";
			// command = "echo -n "+filename+",Dynamic Programming Approach with Modified Kruskal," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./DynamicProgrammingCutsModKrus graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			// system(command.c_str());

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

			// filename = "graph"  +to_string(i) + sample[j] + "half.txt";
			// command = "echo -n "+filename+",Best Cuts with Kruskal Algorithim," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./bestCutsModKrus graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			// system(command.c_str());

			// filename = "graph"  +to_string(i) + sample[j] + "full.txt";
			// command = "echo -n "+filename+",Best Cuts with Kruskal Algorithim," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./bestCutsModKrus graphs/"+filename + " " +to_string(k)+ " >> "+ output;

			// system(command.c_str());


			// filename = "graph"  +to_string(i) + sample[j] + "full.txt";
			// command  = "./graphToMatrix graphs/"+filename;
			// system(command.c_str());

			// filename = "graph"  +to_string(i) + sample[j] + "half.txt";
			// command  = "./graphToMatrix graphs/"+filename;
			// system(command.c_str());


			// filename = "graph"  +to_string(i) + sample[j] + "half.txt";
			// filename2 = "graph"  +to_string(i) + sample[j] + "halfEigen.txt";
			// command = "echo -n "+filename+",Segmantation Clustring with Modified Kruskal Algorithim," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./eigenToGraphs graphs/"+filename + " matrices/" + filename2 + " >> "+ output;

			// system(command.c_str());

			// filename = "graph"  +to_string(i) + sample[j] + "full.txt";
			// filename2 = "graph"  +to_string(i) + sample[j] + "fullEigen.txt";
			// command = "echo -n "+filename+",Segmantation Clustring with Modified Kruskal Algorithim," + to_string(k)+", >> "+output;
			// system(command.c_str());
			// command  = "./eigenToGraphs graphs/"+filename + " matrices/" + filename2+" >> "+ output;

			// system(command.c_str());

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
