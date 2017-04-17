/**
 * \class	Mesh
 * \brief	Stores and manipulates the DCEL decomposition of the affine Grassmannian.
 * \author	Matthew L. Wright
 * \date	March 2014
 */

#ifndef __DCEL_Mesh_H__
#define __DCEL_Mesh_H__

//forward declarations
class Face;
class Halfedge;
class Vertex;

#include "anchor.h"
#include "point.h"

#include <boost/multiprecision/cpp_int.hpp>
typedef boost::multiprecision::cpp_rational exact;

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
typedef boost::property<boost::edge_weight_t, unsigned long> EdgeWeightProperty;
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeWeightProperty > Graph;


#include <vector>
#include <set>


class Mesh
{
    friend class PersistenceUpdater; //allow PersistenceUpdater access to private variables in Mesh

    public:
        Mesh(const std::vector<double>& xg, const std::vector<exact>& xe, const std::vector<double>& yg, const std::vector<exact>& ye, int verbosity);
            //constructor; sets up bounding box (with empty interior) for the affine Grassmannian
		
        ~Mesh();	//destructor
		
        void build_full_arrangement(std::vector<Point>& pts);      //builds the full arrangement (in entire plane)
        void build_half_arrangement(std::vector<Point>& pts);      //builds the half arrangement (in x >=0 half-plane)

        void build_graph(Graph& dual_graph, int &num_edges); //builds the dual graph of the arrangement


        unsigned num_faces();   //returns the number of 2-cells, and thus the number of barcode templates, in the arrangement
		
        //TESTING ONLY
        void print_stats(); //prints a summary of the arrangement information, such as the number of anchors, vertices, halfedges, and faces
		void print();	//prints all the data from the mesh
        void test_consistency();    //attempts to find inconsistencies in the DCEL arrangement
		
        //references to vectors of grade values
        const std::vector<double>& x_grades;   //floating-point values for x-grades
        const std::vector<exact>& x_exact;     //exact values for all x-grades
        const std::vector<double>& y_grades;   //floating-point values for y-grades
        const std::vector<exact>& y_exact;     //exact values for all y-grades


        //these are necessary for comparisons, but should they really be static members of Mesh???
        static double epsilon;
        static bool almost_equal(const double a, const double b);

    private:
      //data structures
        std::vector<Vertex*> vertices;		//all vertices in the mesh
		std::vector<Halfedge*> halfedges;	//all halfedges in the mesh
		std::vector<Face*> faces;		//all faces in the mesh
		
		const double INFTY;

        std::set<Anchor*, Anchor_LeftComparator_Full> all_anchors_full;	//set of Anchors that are represented in the mesh, ordered by position of curve along left side of the arrangement, from bottom to top
        std::set<Anchor*, Anchor_LeftComparator_Half> all_anchors_half;	//set of Anchors that are represented in the mesh, ordered by position of curve along left side of the arrangement, from bottom to top
		
        Halfedge* topleft;			//pointer to Halfedge that points down from top left corner (0,infty)
        Halfedge* topright;         //pointer to Halfedge that points down from the top right corner (infty,infty)
        Halfedge* bottomleft;       //pointer to Halfedge that points up from bottom left corner (0,-infty)
        Halfedge* bottomright;      //pointer to Halfedge that points up from bottom right corner (infty,-infty)
		
        std::vector<Halfedge*> vertical_line_query_list; //stores a pointer to the rightmost Halfedge of the "top" line of each unique slope, ordered from small slopes to big slopes (each Halfedge points to Anchor and Face for vertical-line queries)

        const int verbosity;			//controls display of output, for debugging


      //functions for creating the arrangement
        void build_full_interior();
            //builds the interior of DCEL arrangement using a version of the Bentley-Ottmann algorithm
            //precondition: all achors have been stored via find_anchors()

        void build_half_interior();

        Halfedge* insert_vertex(Halfedge* edge, double x, double y);	//inserts a new vertex on the specified edge, with the specified coordinates, and updates all relevant pointers
        Halfedge* create_edge_left(Halfedge* edge, Anchor *anchor);    //creates the first pair of Halfedges in an anchor line, anchored on the left edge of the strip


      //functions for testing
        unsigned HID(Halfedge* h);		//halfedge ID, for printing and debugging
        unsigned FID(Face* f);		//face ID, for printing and debugging
        unsigned VID(Vertex* v);    //vertex ID, for printing and debugging


      //struct to hold a future intersection event
        struct Crossing {
            Anchor* a;     //pointer to one line
            Anchor* b;     //pointer to the other line -- must ensure that line for anchor a is below line for anchor b just before the crossing point!!!!!
            double x;   //x-coordinate of intersection point (floating-point)
            Mesh* m;    //pointer to the mesh, so the Crossing has access to the vectors x_grades, x_exact, y_grades, and y_exact

            Crossing(Anchor* a, Anchor* b, Mesh* m);  //precondition: Anchors a and b must be comparable
            bool x_equal(const Crossing* other) const;  //returns true iff this Crossing has (exactly) the same x-coordinate as other Crossing
        };

      //comparator class for ordering crossings: first by x (left to right); for a given x, then by y (low to high)
        struct CrossingComparator {
            bool operator()(const Crossing* c1, const Crossing* c2) const;	//returns true if c1 comes after c2
        };

};//end class Mesh

#endif // __DCEL_Mesh_H__

