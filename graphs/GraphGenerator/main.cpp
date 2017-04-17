#include "dcel/mesh.h"
#include "point.h"

#include <qdebug.h>
#include <QFile>
#include <QTextStream>

#include <boost/multiprecision/cpp_int.hpp>
typedef boost::multiprecision::cpp_rational exact;

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
typedef boost::property<boost::edge_weight_t, unsigned long> EdgeWeightProperty;
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeWeightProperty > Graph;


#include <cstdlib>      //for atoi
#include <stdlib.h>     //for srand, rand
#include <time.h>       //for time


int main(int argc, char *argv[])
{
    //parameters
    if(argc < 2)
    {
        qDebug() << "Usage: " << argv[0] << "num_points";
        return 1;
    }
    int num_pts = std::atoi(argv[1]);

    //parameters
    int min_coord = -20;
    int max_coord = 20;
    int max_weight = 1000;
    int repeat = 3;

    //set up the grid
    std::vector<double> x_grades;
    std::vector<double> y_grades;
    std::vector<exact> x_exact;
    std::vector<exact> y_exact;

    for(int i = min_coord; i <= max_coord; i++)
    {
        x_grades.push_back(i);
        y_grades.push_back(i);
        x_exact.push_back(i);
        y_exact.push_back(i);
    }

    //seed the random number generator
    srand(time(NULL));

    //generate 3 sets of graphs
    char ind = 'A';
    for(int s = 0; s < repeat; s++)
    {

        //generate a list of random points
        qDebug() << "Generating" << num_pts << "random points";
        std::vector<Point> pts;
        for(int i = 0; i < num_pts; i++)
        {
            int x = rand() % (max_coord - min_coord + 1);
            int y = rand() % (max_coord - min_coord + 1);
            int w = rand() % max_weight + 1;
            pts.push_back(Point(x, y, w));
            qDebug() << "  (" << x_grades[x] << "," << y_grades[y] << "), weight =" << w;
        }

        //create the full graph
        Mesh* full_arr = new Mesh(x_grades, x_exact, y_grades, y_exact, 2);    //delete later
        full_arr->build_full_arrangement(pts);
        Graph full_graph;
        int full_edges=0;
        full_arr->build_graph(full_graph, full_edges);

        //output the full graph
        qDebug() << "creating the full graph... ";
        typedef boost::graph_traits<Graph>::edge_iterator edge_iterator;
        std::pair<edge_iterator, edge_iterator> ei = boost::edges(full_graph);
        QString filename = QString("graph") + QString::number(num_pts) + QString(ind) + QString("full.txt");
        QFile file(filename);
        if(file.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text))
        {
            QTextStream stream(&file);
            stream << "p   edge " << full_arr->num_faces() << "  " << full_edges << endl;
            for(edge_iterator it = ei.first; it != ei.second; ++it)
                stream << "e  " << boost::source(*it, full_graph) << "  " << boost::target(*it, full_graph) << "  " << boost::get(boost::edge_weight_t(), full_graph, *it) << endl;
        }

        //create the half graph
        Mesh* half_arr = new Mesh(x_grades, x_exact, y_grades, y_exact, 2);    //delete later
        half_arr->build_half_arrangement(pts);
        Graph half_graph;
        int half_edges=0;
        half_arr->build_graph(half_graph, half_edges);

        //output the full graph
        qDebug() << "creating the half graph...";
        typedef boost::graph_traits<Graph>::edge_iterator edge_iterator;
        std::pair<edge_iterator, edge_iterator> ei2 = boost::edges(half_graph);
        QString filename2 = QString("graph") + QString::number(num_pts) + QString(ind) + QString("half.txt");
        QFile file2(filename2);
        if(file2.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text))
        {
            QTextStream stream(&file2);
            stream << "p   edge " << half_arr->num_faces() << "  " << half_edges << endl;
            for(edge_iterator it = ei2.first; it != ei2.second; ++it)
                stream << "e  " << boost::source(*it, half_graph) << "  " << boost::target(*it, half_graph) << "  " << boost::get(boost::edge_weight_t(), half_graph, *it) << endl;
        }

        //clean up
        delete full_arr;
        delete half_arr;

        //increment the index
        ind++;
    }

    //done
    return 0;
}//end main()
