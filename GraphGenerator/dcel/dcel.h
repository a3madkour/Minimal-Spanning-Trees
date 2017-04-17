/**
 * \brief	Vertex, Halfedge, and Face classes for building a DCEL arrangement
 *
 * A DCEL is a doubly connected edge list, see e.g. de Berg et. al.'s Computational
 * Geometry: Algorithms and Applications for background.
 * 
 * \author	Matthew L. Wright
 * \date	March 2014
 */

#ifndef __DCEL_H__
#define __DCEL_H__

//forward declarations
class Vertex;
class Halfedge;
class Face;
class Anchor;

#include "barcode_template.h"

#include <QDebug>


class Vertex
{
    public:
        Vertex(double x_coord, double y_coord);	//constructor, sets (x, y)-coordinates of the vertex

        void set_incident_edge(Halfedge* edge);		//set the incident edge
        Halfedge* get_incident_edge();			//get the incident edge

        double get_x();		//get the x-coordinate
        double get_y();		//get the y-coordinate

        friend QDebug& operator<<(QDebug& qd, const Vertex& v);	//for printing the vertex

    private:
        Halfedge* incident_edge;	//pointer to one edge incident to this vertex
        double x;			//x-coordinate of this vertex
        double y;			//y-coordinate of this vertex

};//end class Vertex


class Halfedge
{
    public:
        Halfedge(Vertex* v, Anchor* p); //constructor, requires origin vertex as well as Anchor corresponding to this halfedge (Anchor never changes)
        Halfedge();                     //constructor for a null Halfedge

        void set_twin(Halfedge* e); //set the twin halfedge
        Halfedge* get_twin() const; //get the twin halfedge

        void set_next(Halfedge* e); //set the next halfedge in the boundary of the face that this halfedge borders
        Halfedge* get_next() const; //get the next halfedge

        void set_prev(Halfedge* e); //set the previous halfedge in the boundary of the face that this halfedge borders
        Halfedge* get_prev() const; //get the previous halfedge

        void set_origin(Vertex* v); //set the origin vertex
        Vertex* get_origin() const; //get the origin vertex

        void set_face(Face* f);     //set the face that this halfedge borders
        Face* get_face() const;     //get the face that this halfedge borders

        Anchor* get_anchor() const; //get the Anchor

        friend QDebug& operator<<(QDebug& qd, const Halfedge& e);	//for printing the halfedge

    private:
        Vertex* origin;		//pointer to the vertex from which this halfedge originates
        Halfedge* twin;		//pointer to the halfedge that, together with this halfedge, make one edge
        Halfedge* next;		//pointer to the next halfedge around the boundary of the face to the right of this halfedge
        Halfedge* prev;		//pointer to the previous halfedge around the boundary of the face to the right of this halfedge
        Face* face;		    //pointer to the face to the right of this halfedge
        Anchor* anchor;		//stores the coordinates of the anchor corresponding to this halfedge

};//end class Halfedge


class Face
{
    public:
        Face(Halfedge* e);      //constructor: requires pointer to a boundary halfedge
        ~Face();                //destructor: destroys barcode template

        void set_boundary(Halfedge* e);     //set the pointer to a halfedge on the boundary of this face
        Halfedge* get_boundary();           //get the (pointer to the) boundary halfedge

        BarcodeTemplate& get_barcode();         //returns a reference to the barcode template stored in this cell
        void set_barcode(BarcodeTemplate& bt);  //stores (a copy of) the specified barcode template in this cell

        bool has_been_visited();            //true iff cell has been visited in the vineyard-update process (so that we can distinguish a cell with an empty barcode from an unvisited cell)
        void mark_as_visited();             //marks this cell as visited

        friend QDebug& operator<<(QDebug& os, const Face& f);	//for printing the face

    private:
        Halfedge* boundary;     //pointer to one halfedge in the boundary of this cell
        BarcodeTemplate dbc;    //barcode template stored in this cell
        bool visited;           //initially false, set to true after this cell has been visited in the vineyard-update process (so that we can distinguish a cell with an empty barcode from an unvisited cell)

};//end class Face


#endif // __DCEL_H__
