/* implementation of Mesh class
 * Stores and manipulates the DCEL decomposition of the affine Grassmannian.
 */

#include "mesh.h"

#include "dcel.h"
#include "../point.h"

#include <QDebug>
#include <QTime>



#include <limits>	//necessary for infinity
#include <queue>    //for std::priority_queue


// Mesh constructor; sets up bounding box (with empty interior) for the affine Grassmannian
Mesh::Mesh(const std::vector<double> &xg, const std::vector<exact> &xe, const std::vector<double> &yg, const std::vector<exact> &ye, int verbosity) :
    x_grades(xg), x_exact(xe), y_grades(yg), y_exact(ye),
    INFTY(std::numeric_limits<double>::infinity()),
    verbosity(verbosity)
{
    //create vertices
    vertices.push_back( new Vertex(-INFTY, INFTY) );         //index 0
    vertices.push_back( new Vertex(INFTY, INFTY) );     //index 1
    vertices.push_back( new Vertex(INFTY, -INFTY) );    //index 2
    vertices.push_back( new Vertex(-INFTY, -INFTY) );        //index 3

    //create halfedges
    for(int i=0; i<4; i++)
    {
        halfedges.push_back( new Halfedge( vertices[i], NULL) );		//index 0, 2, 4, 6 (inside halfedges)
        halfedges.push_back( new Halfedge( vertices[(i+1)%4], NULL) );		//index 1, 3, 5, 7 (outside halfedges)
        halfedges[2*i]->set_twin( halfedges[2*i+1] );
        halfedges[2*i+1]->set_twin( halfedges[2*i] );
    }

    topleft = halfedges[7];     //remember this halfedge to make curve insertion easier
    topright = halfedges[2];    //remember this halfedge for starting the path that we use to find edge weights
    bottomleft = halfedges[6];  //remember these halfedges
    bottomright = halfedges[3]; //    for the Bentley-Ottmann algorithm

    //set pointers on vertices
    for(int i=0; i<4; i++)
    {
        vertices[i]->set_incident_edge( halfedges[2*i] );
    }

    //create face
    faces.push_back( new Face( halfedges[0] ) );

    //set the remaining pointers on the halfedges
    for(int i=0; i<4; i++)
    {
        Halfedge* inside = halfedges[2*i];
        inside->set_next( halfedges[(2*i+2)%8] );
        inside->set_prev( halfedges[(2*i+6)%8] );
        inside->set_face( faces[0] );

        Halfedge* outside = halfedges[2*i+1];
        outside->set_next( halfedges[(2*i+7)%8] );
        outside->set_prev( halfedges[(2*i+3)%8] );
    }
}//end constructor

//destructor
Mesh::~Mesh()
{
    for(std::vector<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); ++it)
        delete (*it);

    for(std::vector<Halfedge*>::iterator it = halfedges.begin(); it != halfedges.end(); ++it)
        delete (*it);

    for(std::vector<Face*>::iterator it = faces.begin(); it != faces.end(); ++it)
        delete (*it);

    for(std::set<Anchor*>::iterator it = all_anchors_full.begin(); it != all_anchors_full.end(); ++it)
        delete (*it);

    for(std::set<Anchor*>::iterator it = all_anchors_half.begin(); it != all_anchors_half.end(); ++it)
        delete (*it);
}//end destructor

//precondition: the constructor has already created the boundary of the arrangement
void Mesh::build_full_arrangement(std::vector<Point>& pts)
{
    //FOR GRAPH TESTS, WE REGARD ALL GIVEN POINTS AS ANCHORS!
    for(std::vector<Point>::iterator it = pts.begin(); it != pts.end(); ++it)
    {
//        qDebug() << "adding point: " << x_grades[it->x] << y_grades[it->y];
        Anchor* anch = new Anchor(it->x, it->y);
        anch->set_weight(it->w);
        all_anchors_full.insert(anch);
    }

    //now that we have all the anchors, we can build the interior of the arrangement
    build_full_interior();
//    print_stats();
}

//precondition: the constructor has already created the boundary of the arrangement
void Mesh::build_half_arrangement(std::vector<Point>& pts)
{
    //FOR GRAPH TESTS, WE REGARD ALL GIVEN POINTS AS ANCHORS!
    for(std::vector<Point>::iterator it = pts.begin(); it != pts.end(); ++it)
    {
//        qDebug() << "adding point: " << x_grades[it->x] << y_grades[it->y];
        Anchor* anch = new Anchor(it->x, it->y);
        anch->set_weight(it->w);
        all_anchors_half.insert(anch);
    }

    //now that we have all the anchors, we can build the interior of the arrangement
    build_half_interior();
//    print_stats();
}

//builds the dual graph of the arrangement
//precondition: graph is empty
void Mesh::build_graph(Graph& dual_graph, int& num_edges)
{
    //make a map for reverse-lookup of face indexes by face pointers -- Is this really necessary? Can I avoid doing this?
    std::map<Face*, unsigned> face_indexes;
    for(unsigned i=0; i<faces.size(); i++)
        face_indexes.insert( std::pair<Face*, unsigned>(faces[i], i));

    //loop over all faces
    for(unsigned i=0; i<faces.size(); i++)
    {
        //consider all neighbors of this faces
        Halfedge* boundary = (faces[i])->get_boundary();
        Halfedge* current = boundary;
        do
        {
            //find index of neighbor
            Face* neighbor = current->get_twin()->get_face();
            if(neighbor != NULL)
            {
                std::map<Face*, unsigned>::iterator it = face_indexes.find(neighbor);
                unsigned j = it->second;

                //if i < j, then create an (undirected) edge between these faces
                if(i < j)
                {
                    boost::add_edge(i, j, current->get_anchor()->get_weight(), dual_graph);
                    num_edges++;
                }
            }
              //move to the next neighbor
            current = current->get_next();
        }while(current != boundary);
    }

}//end build_graph()


//function to build the arrangement using a version of the Bentley-Ottmann algorithm, given all Anchors
//preconditions:
//   all Anchors are in a list, ordered by Anchor_LeftComparator
//   boundary of the mesh is created (as in the mesh constructor)
void Mesh::build_full_interior()
{
    if(verbosity >= 4)
    {
        QDebug qd = qDebug().nospace();
        qd << "BUILDING FULL ARRANGEMENT:  Points sorted for left edge: ";
        for(std::set<Anchor*, Anchor_LeftComparator_Full>::iterator it = all_anchors_full.begin(); it != all_anchors_full.end(); ++it)
            qd << "(" << x_grades[(*it)->get_x()] << "," << y_grades[(*it)->get_y()] << ") ";
    }

  // DATA STRUCTURES

    //data structure for ordered list of lines
    std::vector<Halfedge*> lines;
    lines.reserve(all_anchors_full.size());

    //data structure for queue of future intersections
    std::priority_queue< Crossing*, std::vector<Crossing*>, CrossingComparator > crossings;

    //data structure for all pairs of Anchors whose potential crossings have been considered
    typedef std::pair<Anchor*,Anchor*> Anchor_pair;
    std::set< Anchor_pair > considered_pairs;

  // PART 1: INSERT VERTICES AND EDGES ALONG LEFT EDGE OF THE ARRANGEMENT
    if(verbosity >= 6) { qDebug() << "PART 1: LEFT EDGE OF ARRANGEMENT"; }

    //for each Anchor, create vertex and associated halfedges, anchored on the left edge of the strip
    Halfedge* leftedge = bottomleft;
    for(std::set<Anchor*, Anchor_LeftComparator_Full>::iterator it = all_anchors_full.begin(); it != all_anchors_full.end(); ++it)
    {
        Anchor* cur_anchor = *it;

        if(verbosity >= 6) { qDebug() << "  Processing Anchor" << cur_anchor << "at (" << x_grades[cur_anchor->get_x()] << "," << y_grades[cur_anchor->get_y()] << ")"; }

        leftedge = insert_vertex(leftedge, -INFTY, 0);    //set leftedge to edge that will follow the new edge

        //now insert new edge at origin vertex of leftedge
        Halfedge* new_edge = create_edge_left(leftedge, cur_anchor);

        //remember Halfedge corresponding to this Anchor
        lines.push_back(new_edge);

        //remember relative position of this Anchor
        cur_anchor->set_position(lines.size() - 1);

        //remember line associated with this Anchor
        cur_anchor->set_line(new_edge);
    }

    //for each pair of consecutive lines, if they intersect, store the intersection
    for(unsigned i = 0; i+1 < lines.size(); i++)
    {
        Anchor* a = lines[i]->get_anchor();
        Anchor* b = lines[i+1]->get_anchor();
        if( a->get_x() != b->get_x() )    //then the lines are not parallel
            crossings.push(new Crossing(a, b, this));

        //remember that we have now considered this intersection
        considered_pairs.insert(Anchor_pair(a,b));
    }


  // PART 2: PROCESS INTERIOR INTERSECTIONS
    //    order: x left to right; for a given x, then y low to high
    if(verbosity >= 6) { qDebug() << "PART 2: PROCESSING INTERIOR INTERSECTIONS\n"; }

    int status_counter = 0;
    int status_interval = 10000;    //controls frequency of output

    //current position of sweep line
    Crossing* sweep = NULL;

    while(!crossings.empty())
    {
        //get the next intersection from the queue
        Crossing* cur = crossings.top();
        crossings.pop();

        //process the intersection
        sweep = cur;
        unsigned first_pos = cur->a->get_position();   //most recent edge in the curve corresponding to Anchor a
        unsigned last_pos = cur->b->get_position();   //most recent edge in the curve corresponding to Anchor b

        if(verbosity >= 6) { qDebug() << " next intersection: Anchor" << cur->a << " (pos" << first_pos << "), Anchor" << cur->b << " (pos" << last_pos << ") at x =" << cur->x; }

        if(last_pos != first_pos + 1)
        {
            qDebug() << "ERROR: intersection between non-consecutive curves [1]: x = " << sweep->x << "\n";
            throw std::exception();
        }

        //find out if more than two curves intersect at this point
        while( !crossings.empty() && sweep->x_equal(crossings.top()) && (cur->b == crossings.top()->a) )
        {
            cur = crossings.top();
            crossings.pop();

            if(cur->b->get_position() != last_pos + 1)
            {
                qDebug() << "ERROR: intersection between non-consecutive curves [2]\n";
                throw std::exception();
            }

            last_pos++; //last_pos = cur->b->get_position();

            if(verbosity >= 6) { qDebug() << " |---also intersects Anchor" << cur->b << "(" << last_pos << ")"; }
        }

        //compute y-coordinate of intersection
        double intersect_y = x_grades[sweep->a->get_x()]*(sweep->x) - y_grades[sweep->a->get_y()];

        if(verbosity >= 8) { qDebug() << "  found intersection between" << (last_pos - first_pos + 1) << "edges at x =" << sweep->x << ", y =" << intersect_y; }

        //create new vertex
        Vertex* new_vertex = new Vertex(sweep->x, intersect_y);
        vertices.push_back(new_vertex);

        //anchor edges to vertex and create new face(s) and edges	//TODO: check this!!!
        Halfedge* prev_new_edge = NULL;                 //necessary to remember the previous new edge at each interation of the loop
        Halfedge* first_incoming = lines[first_pos];   //necessary to remember the first incoming edge
        Halfedge* prev_incoming = NULL;                 //necessary to remember the previous incoming edge at each iteration of the loop
        for(unsigned cur_pos = first_pos; cur_pos <= last_pos; cur_pos++)
        {
            //anchor edge to vertex
            Halfedge* incoming = lines[cur_pos];
            incoming->get_twin()->set_origin(new_vertex);

            //create next pair of twin halfedges along the current curve (i.e. curves[incident_edges[i]] )
            Halfedge* new_edge = new Halfedge(new_vertex, incoming->get_anchor());	//points AWAY FROM new_vertex
            halfedges.push_back(new_edge);
            Halfedge* new_twin = new Halfedge(NULL, incoming->get_anchor());		//points TOWARDS new_vertex
            halfedges.push_back(new_twin);

            //update halfedge pointers
            new_edge->set_twin(new_twin);
            new_twin->set_twin(new_edge);

            if(cur_pos == first_pos)    //then this is the first iteration of the loop
            {
                new_twin->set_next( lines[last_pos]->get_twin() );
                lines[last_pos]->get_twin()->set_prev(new_twin);

                new_twin->set_face( lines[last_pos]->get_twin()->get_face() );
            }
            else    //then this is not the first iteration of the loop, so close up a face and create a new face
            {
                incoming->set_next( prev_incoming->get_twin() );
                incoming->get_next()->set_prev(incoming);

                Face* new_face = new Face(new_twin);
                faces.push_back(new_face);

                new_twin->set_face(new_face);
                prev_new_edge->set_face(new_face);

                new_twin->set_next(prev_new_edge);
                prev_new_edge->set_prev(new_twin);
            }

            //remember important halfedges for the next iteration of the loop
            prev_incoming = incoming;
            prev_new_edge = new_edge;

            if(cur_pos == last_pos)  //then this is the last iteration of loop
            {
                new_edge->set_prev(first_incoming);
                first_incoming->set_next(new_edge);

                new_edge->set_face( first_incoming->get_face() );
            }

            //update lines vector
            lines[cur_pos] = new_edge; //the portion of this vector [first_pos, last_pos] must be reversed after this loop is finished!

            //remember position of this Anchor
            new_edge->get_anchor()->set_position(last_pos - (cur_pos - first_pos));
        }

        //update lines vector: flip portion of vector [first_pos, last_pos]
        for(unsigned i = 0; i < (last_pos - first_pos + 1)/2; i++)
        {
            //swap curves[first_pos + i] and curves[last_pos - i]
            Halfedge* temp = lines[first_pos + i];
            lines[first_pos + i] = lines[last_pos - i];
            lines[last_pos - i] = temp;
        }

        //find new intersections and add them to intersections queue
        if(first_pos > 0)   //then consider lower intersection
        {
            Anchor* a = lines[first_pos-1]->get_anchor();
            Anchor* b = lines[first_pos]->get_anchor();

            if(considered_pairs.find(Anchor_pair(a,b)) == considered_pairs.end()
                    && considered_pairs.find(Anchor_pair(b,a)) == considered_pairs.end() )	//then this pair has not yet been considered
            {
                considered_pairs.insert(Anchor_pair(a,b));
                    crossings.push(new Crossing(a, b, this));
            }
        }

        if(last_pos + 1 < lines.size())    //then consider upper intersection
        {
            Anchor* a = lines[last_pos]->get_anchor();
            Anchor* b = lines[last_pos+1]->get_anchor();

            if( considered_pairs.find(Anchor_pair(a,b)) == considered_pairs.end()
                    && considered_pairs.find(Anchor_pair(b,a)) == considered_pairs.end() )	//then this pair has not yet been considered
            {
                considered_pairs.insert(Anchor_pair(a,b));
                    crossings.push(new Crossing(a, b, this));
            }
        }

        //output status
        if(verbosity >= 4)
        {
            status_counter++;
            if(status_counter % status_interval == 0)
                qDebug() << "      processed" << status_counter << "intersections; sweep position =" << sweep;
        }
    }//end while

  // PART 3: INSERT VERTICES ON RIGHT EDGE OF ARRANGEMENT AND CONNECT EDGES
    if(verbosity >= 6) { qDebug() << "PART 3: RIGHT EDGE OF THE ARRANGEMENT"; }

    Halfedge* rightedge = bottomright; //need a reference halfedge along the right side of the strip
    unsigned cur_x = 0;      //keep track of discrete x-coordinate of last Anchor whose line was connected to right edge (x-coordinate of Anchor is slope of line)

    //connect each line to the right edge of the arrangement (at x = INFTY)
    //    requires creating a vertex for each unique slope (i.e. Anchor x-coordinate)
    //    lines that have the same slope m are "tied together" at the same vertex, with coordinates (INFTY, Y)
    //    where Y = INFTY if m is positive, Y = -INFTY if m is negative, and Y = 0 if m is zero
    for(unsigned cur_pos = 0; cur_pos < lines.size(); cur_pos++)
    {
        Halfedge* incoming = lines[cur_pos];
        Anchor* cur_anchor = incoming->get_anchor();

        if(cur_anchor->get_x() > cur_x || cur_pos == 0)    //then create a new vertex for this line
        {
            cur_x = cur_anchor->get_x();

            double Y = INFTY;               //default, for lines with positive slope
            if(x_grades[cur_x] < 0)
                Y = -1*Y;                   //for lines with negative slope
            else if(x_grades[cur_x] == 0)
                Y = 0;                      //for horizontal lines

            rightedge = insert_vertex( rightedge, INFTY, Y );
        }
        else    //no new vertex required, but update previous entry for vertical-line queries
            vertical_line_query_list.pop_back();

        //store Halfedge for vertical-line queries
        vertical_line_query_list.push_back(incoming->get_twin());

        //connect current line to the most-recently-inserted vertex
        Vertex* cur_vertex = rightedge->get_origin();
        incoming->get_twin()->set_origin(cur_vertex);

        //update halfedge pointers
        incoming->set_next(rightedge->get_twin()->get_next());
        incoming->get_next()->set_prev(incoming);

        incoming->get_next()->set_face(incoming->get_face());   //only necessary if incoming->get_next() is along the right side of the strip

        incoming->get_twin()->set_prev(rightedge->get_twin());
        rightedge->get_twin()->set_next(incoming->get_twin());

        rightedge->get_twin()->set_face(incoming->get_twin()->get_face());
    }
}//end build_full_interior()

//function to build the arrangement using a version of the Bentley-Ottmann algorithm, given all Anchors
//preconditions:
//   all Anchors are in a list, ordered by Anchor_LeftComparator
//   boundary of the mesh is created (as in the mesh constructor)
void Mesh::build_half_interior()
{
    if(verbosity >= 4)
    {
        QDebug qd = qDebug().nospace();
        qd << "BUILDING HALF ARRANGEMENT:  Points sorted for left edge: ";
        for(std::set<Anchor*, Anchor_LeftComparator_Half>::iterator it = all_anchors_half.begin(); it != all_anchors_half.end(); ++it)
            qd << "(" << x_grades[(*it)->get_x()] << "," << y_grades[(*it)->get_y()] << ") ";
    }

  // DATA STRUCTURES

    //data structure for ordered list of lines
    std::vector<Halfedge*> lines;
    lines.reserve(all_anchors_half.size());

    //data structure for queue of future intersections
    std::priority_queue< Crossing*, std::vector<Crossing*>, CrossingComparator > crossings;

    //data structure for all pairs of Anchors whose potential crossings have been considered
    typedef std::pair<Anchor*,Anchor*> Anchor_pair;
    std::set< Anchor_pair > considered_pairs;

  // PART 1: INSERT VERTICES AND EDGES ALONG LEFT EDGE OF THE ARRANGEMENT
    if(verbosity >= 6) { qDebug() << "PART 1: LEFT EDGE OF ARRANGEMENT"; }

    //for each Anchor, create vertex and associated halfedges, anchored on the left edge of the strip
    Halfedge* leftedge = bottomleft;
    unsigned prev_y = std::numeric_limits<unsigned>::max();
    for(std::set<Anchor*, Anchor_LeftComparator_Half>::iterator it = all_anchors_half.begin(); it != all_anchors_half.end(); ++it)
    {
        Anchor* cur_anchor = *it;

        if(verbosity >= 6) { qDebug() << "  Processing Anchor" << cur_anchor << "at (" << x_grades[cur_anchor->get_x()] << "," << y_grades[cur_anchor->get_y()] << ")"; }

        if(cur_anchor->get_y() != prev_y)	//then create new vertex
        {
            double dual_point_y_coord = -1*y_grades[cur_anchor->get_y()];  //point-line duality requires multiplying by -1
            leftedge = insert_vertex(leftedge, 0, dual_point_y_coord);    //set leftedge to edge that will follow the new edge
            prev_y = cur_anchor->get_y();  //remember the discrete y-index
        }

        //now insert new edge at origin vertex of leftedge
        Halfedge* new_edge = create_edge_left(leftedge, cur_anchor);

        //remember Halfedge corresponding to this Anchor
        lines.push_back(new_edge);

        //remember relative position of this Anchor
        cur_anchor->set_position(lines.size() - 1);

        //remember line associated with this Anchor
        cur_anchor->set_line(new_edge);
    }

    //for each pair of consecutive lines, if they intersect, store the intersection
    for(unsigned i = 0; i+1 < lines.size(); i++)
    {
        Anchor* a = lines[i]->get_anchor();
        Anchor* b = lines[i+1]->get_anchor();
        if( a->comparable(b) )    //then the lines are not parallel
            crossings.push(new Crossing(a, b, this));

        //remember that we have now considered this intersection
        considered_pairs.insert(Anchor_pair(a,b));
    }


  // PART 2: PROCESS INTERIOR INTERSECTIONS
    //    order: x left to right; for a given x, then y low to high
    if(verbosity >= 6) { qDebug() << "PART 2: PROCESSING INTERIOR INTERSECTIONS\n"; }

    int status_counter = 0;
    int status_interval = 10000;    //controls frequency of output

    //current position of sweep line
    Crossing* sweep = NULL;

    while(!crossings.empty())
    {
        //get the next intersection from the queue
        Crossing* cur = crossings.top();
        crossings.pop();

        //process the intersection
        sweep = cur;
        unsigned first_pos = cur->a->get_position();   //most recent edge in the curve corresponding to Anchor a
        unsigned last_pos = cur->b->get_position();   //most recent edge in the curve corresponding to Anchor b

        if(verbosity >= 6) { qDebug() << " next intersection: Anchor" << cur->a << " (pos" << first_pos << "), Anchor" << cur->b << " (pos" << last_pos << ") at x =" << cur->x; }

        if(last_pos != first_pos + 1)
        {
            qDebug() << "ERROR: intersection between non-consecutive curves [1]: x = " << sweep->x << "\n";
            throw std::exception();
        }

        //find out if more than two curves intersect at this point
        while( !crossings.empty() && sweep->x_equal(crossings.top()) && (cur->b == crossings.top()->a) )
        {
            cur = crossings.top();
            crossings.pop();

            if(cur->b->get_position() != last_pos + 1)
            {
                qDebug() << "ERROR: intersection between non-consecutive curves [2]\n";
                throw std::exception();
            }

            last_pos++; //last_pos = cur->b->get_position();

            if(verbosity >= 6) { qDebug() << " |---also intersects Anchor" << cur->b << "(" << last_pos << ")"; }
        }

        //compute y-coordinate of intersection
        double intersect_y = x_grades[sweep->a->get_x()]*(sweep->x) - y_grades[sweep->a->get_y()];

        if(verbosity >= 8) { qDebug() << "  found intersection between" << (last_pos - first_pos + 1) << "edges at x =" << sweep->x << ", y =" << intersect_y; }

        //create new vertex
        Vertex* new_vertex = new Vertex(sweep->x, intersect_y);
        vertices.push_back(new_vertex);

        //anchor edges to vertex and create new face(s) and edges	//TODO: check this!!!
        Halfedge* prev_new_edge = NULL;                 //necessary to remember the previous new edge at each interation of the loop
        Halfedge* first_incoming = lines[first_pos];   //necessary to remember the first incoming edge
        Halfedge* prev_incoming = NULL;                 //necessary to remember the previous incoming edge at each iteration of the loop
        for(unsigned cur_pos = first_pos; cur_pos <= last_pos; cur_pos++)
        {
            //anchor edge to vertex
            Halfedge* incoming = lines[cur_pos];
            incoming->get_twin()->set_origin(new_vertex);

            //create next pair of twin halfedges along the current curve (i.e. curves[incident_edges[i]] )
            Halfedge* new_edge = new Halfedge(new_vertex, incoming->get_anchor());	//points AWAY FROM new_vertex
            halfedges.push_back(new_edge);
            Halfedge* new_twin = new Halfedge(NULL, incoming->get_anchor());		//points TOWARDS new_vertex
            halfedges.push_back(new_twin);

            //update halfedge pointers
            new_edge->set_twin(new_twin);
            new_twin->set_twin(new_edge);

            if(cur_pos == first_pos)    //then this is the first iteration of the loop
            {
                new_twin->set_next( lines[last_pos]->get_twin() );
                lines[last_pos]->get_twin()->set_prev(new_twin);

                new_twin->set_face( lines[last_pos]->get_twin()->get_face() );
            }
            else    //then this is not the first iteration of the loop, so close up a face and create a new face
            {
                incoming->set_next( prev_incoming->get_twin() );
                incoming->get_next()->set_prev(incoming);

                Face* new_face = new Face(new_twin);
                faces.push_back(new_face);

                new_twin->set_face(new_face);
                prev_new_edge->set_face(new_face);

                new_twin->set_next(prev_new_edge);
                prev_new_edge->set_prev(new_twin);
            }

            //remember important halfedges for the next iteration of the loop
            prev_incoming = incoming;
            prev_new_edge = new_edge;

            if(cur_pos == last_pos)  //then this is the last iteration of loop
            {
                new_edge->set_prev(first_incoming);
                first_incoming->set_next(new_edge);

                new_edge->set_face( first_incoming->get_face() );
            }

            //update lines vector
            lines[cur_pos] = new_edge; //the portion of this vector [first_pos, last_pos] must be reversed after this loop is finished!

            //remember position of this Anchor
            new_edge->get_anchor()->set_position(last_pos - (cur_pos - first_pos));
        }

        //update lines vector: flip portion of vector [first_pos, last_pos]
        for(unsigned i = 0; i < (last_pos - first_pos + 1)/2; i++)
        {
            //swap curves[first_pos + i] and curves[last_pos - i]
            Halfedge* temp = lines[first_pos + i];
            lines[first_pos + i] = lines[last_pos - i];
            lines[last_pos - i] = temp;
        }

        //find new intersections and add them to intersections queue
        if(first_pos > 0)   //then consider lower intersection
        {
            Anchor* a = lines[first_pos-1]->get_anchor();
            Anchor* b = lines[first_pos]->get_anchor();

            if(considered_pairs.find(Anchor_pair(a,b)) == considered_pairs.end()
                    && considered_pairs.find(Anchor_pair(b,a)) == considered_pairs.end() )	//then this pair has not yet been considered
            {
                considered_pairs.insert(Anchor_pair(a,b));
                if( a->comparable(b) )
                    crossings.push(new Crossing(a, b, this));
            }
        }

        if(last_pos + 1 < lines.size())    //then consider upper intersection
        {
            Anchor* a = lines[last_pos]->get_anchor();
            Anchor* b = lines[last_pos+1]->get_anchor();

            if( considered_pairs.find(Anchor_pair(a,b)) == considered_pairs.end()
                    && considered_pairs.find(Anchor_pair(b,a)) == considered_pairs.end() )	//then this pair has not yet been considered
            {
                considered_pairs.insert(Anchor_pair(a,b));
                if(a->comparable(b))
                    crossings.push(new Crossing(a, b, this));
            }
        }

        //output status
        if(verbosity >= 4)
        {
            status_counter++;
            if(status_counter % status_interval == 0)
                qDebug() << "      processed" << status_counter << "intersections; sweep position =" << sweep;
        }
    }//end while

  // PART 3: INSERT VERTICES ON RIGHT EDGE OF ARRANGEMENT AND CONNECT EDGES
    if(verbosity >= 6) { qDebug() << "PART 3: RIGHT EDGE OF THE ARRANGEMENT"; }

    Halfedge* rightedge = bottomright; //need a reference halfedge along the right side of the strip
    unsigned cur_x = 0;      //keep track of discrete x-coordinate of last Anchor whose line was connected to right edge (x-coordinate of Anchor is slope of line)

    //connect each line to the right edge of the arrangement (at x = INFTY)
    //    requires creating a vertex for each unique slope (i.e. Anchor x-coordinate)
    //    lines that have the same slope m are "tied together" at the same vertex, with coordinates (INFTY, Y)
    //    where Y = INFTY if m is positive, Y = -INFTY if m is negative, and Y = 0 if m is zero
    for(unsigned cur_pos = 0; cur_pos < lines.size(); cur_pos++)
    {
        Halfedge* incoming = lines[cur_pos];
        Anchor* cur_anchor = incoming->get_anchor();

        if(cur_anchor->get_x() > cur_x || cur_pos == 0)    //then create a new vertex for this line
        {
            cur_x = cur_anchor->get_x();

            double Y = INFTY;               //default, for lines with positive slope
            if(x_grades[cur_x] < 0)
                Y = -1*Y;                   //for lines with negative slope
            else if(x_grades[cur_x] == 0)
                Y = 0;                      //for horizontal lines

            rightedge = insert_vertex( rightedge, INFTY, Y );
        }
        else    //no new vertex required, but update previous entry for vertical-line queries
            vertical_line_query_list.pop_back();

        //store Halfedge for vertical-line queries
        vertical_line_query_list.push_back(incoming->get_twin());

        //connect current line to the most-recently-inserted vertex
        Vertex* cur_vertex = rightedge->get_origin();
        incoming->get_twin()->set_origin(cur_vertex);

        //update halfedge pointers
        incoming->set_next(rightedge->get_twin()->get_next());
        incoming->get_next()->set_prev(incoming);

        incoming->get_next()->set_face(incoming->get_face());   //only necessary if incoming->get_next() is along the right side of the strip

        incoming->get_twin()->set_prev(rightedge->get_twin());
        rightedge->get_twin()->set_next(incoming->get_twin());

        rightedge->get_twin()->set_face(incoming->get_twin()->get_face());
    }
}//end build_half_interior()


//inserts a new vertex on the specified edge, with the specified coordinates, and updates all relevant pointers
//  i.e. new vertex is between initial and termainal points of the specified edge
//returns pointer to a new halfedge, whose initial point is the new vertex, and that follows the specified edge around its face
Halfedge* Mesh::insert_vertex(Halfedge* edge, double x, double y)
{
	//create new vertex
    Vertex* new_vertex = new Vertex(x, y);
	vertices.push_back(new_vertex);
	
    //get twin and Anchor of this edge
	Halfedge* twin = edge->get_twin();
    Anchor* anchor = edge->get_anchor();
	
	//create new halfedges
    Halfedge* up = new Halfedge(new_vertex, anchor);
	halfedges.push_back(up);
    Halfedge* dn = new Halfedge(new_vertex, anchor);
	halfedges.push_back(dn);
		
	//update pointers
	up->set_next(edge->get_next());
	up->set_prev(edge);
	up->set_twin(twin);
    up->set_face(edge->get_face());
	
	up->get_next()->set_prev(up);
	
	edge->set_next(up);
	edge->set_twin(dn);
	
	dn->set_next(twin->get_next());
	dn->set_prev(twin);
	dn->set_twin(edge);
    dn->set_face(twin->get_face());

	dn->get_next()->set_prev(dn);
	
	twin->set_next(dn);
	twin->set_twin(up);
	
	new_vertex->set_incident_edge(up);
	
	//return pointer to up
	return up;
}//end insert_vertex()

//creates the first pair of Halfedges in an Anchor line, anchored on the left edge of the strip at origin of specified edge
//  also creates a new face (the face below the new edge)
//  CAUTION: leaves NULL: new_edge.next and new_twin.prev
Halfedge* Mesh::create_edge_left(Halfedge* edge, Anchor* anchor)
{
    //create new halfedges
    Halfedge* new_edge = new Halfedge(edge->get_origin(), anchor); //points AWAY FROM left edge
    halfedges.push_back(new_edge);
    Halfedge* new_twin = new Halfedge(NULL, anchor);   //points TOWARDS left edge
    halfedges.push_back(new_twin);

    //create new face
    Face* new_face = new Face(new_edge);
    faces.push_back(new_face);

    //update Halfedge pointers
    new_edge->set_prev(edge->get_prev());
    new_edge->set_twin(new_twin);
    new_edge->set_face(new_face);

    edge->get_prev()->set_next(new_edge);
    edge->get_prev()->set_face(new_face);
    if(edge->get_prev()->get_prev() != NULL)
        edge->get_prev()->get_prev()->set_face(new_face);

    new_twin->set_next(edge);
    new_twin->set_twin(new_edge);
    new_twin->set_face(edge->get_face());

    edge->set_prev(new_twin);

    //return pointer to new_edge
    return new_edge;
}//end create_edge_left()



//returns the number of 2-cells, and thus the number of barcode templates, in the arrangement
unsigned Mesh::num_faces()
{
    return faces.size();
}


//prints a summary of the arrangement information, such as the number of anchors, vertices, halfedges, and faces
void Mesh::print_stats()
{
    qDebug() << "The arrangement contains:" << vertices.size() << "vertices," << halfedges.size() << "halfedges, and" << faces.size() << "faces";
}

//print all the data from the mesh
void Mesh::print()
{
    qDebug() << "  Vertices";
    for(unsigned i=0; i<vertices.size(); i++)
	{
        qDebug() << "    vertex " << i << ": " << *vertices[i] << "; incident edge: " << HID(vertices[i]->get_incident_edge());
	}
	
    qDebug() << "  Halfedges";
    for(unsigned i=0; i<halfedges.size(); i++)
	{
		Halfedge* e = halfedges[i];
		Halfedge* t = e->get_twin();
        qDebug() << "    halfedge " << i << ": " << *(e->get_origin()) << "--" << *(t->get_origin()) << "; ";
        if(e->get_anchor() == NULL)
            qDebug() << "Anchor null; ";
		else
            qDebug() << "Anchor coords (" << e->get_anchor()->get_x() << ", " << e->get_anchor()->get_y() << "); ";
        qDebug() << "twin: " << HID(t) << "; next: " << HID(e->get_next()) << "; prev: " << HID(e->get_prev()) << "; face: " << FID(e->get_face());
	}
	
    qDebug() << "  Faces";
    for(unsigned i=0; i<faces.size(); i++)
	{
        qDebug() << "    face " << i << ": " << *faces[i];
	}
	
/*    qDebug() << "  Outside (unbounded) region: ";
	Halfedge* start = halfedges[1];
	Halfedge* curr = start;
	do{
        qDebug() << *(curr->get_origin()) << "--";
		curr = curr->get_next();
	}while(curr != start);
    qDebug() << "cycle";
*/
//    qDebug() << "  Anchor set: ";
//    std::set<Anchor*>::iterator it;
//    for(it = all_anchors.begin(); it != all_anchors.end(); ++it)
//	{
//        Anchor cur = **it;
//        qDebug() << "(" << cur.get_x() << ", " << cur.get_y() << ") halfedge " << HID(cur.get_line()) << "; ";
//	}
}//end print()

/********** functions for testing **********/

//look up halfedge ID, used in print() for debugging
// HID = halfedge ID
unsigned Mesh::HID(Halfedge* h)
{
    for(unsigned i=0; i<halfedges.size(); i++)
	{
		if(halfedges[i] == h)
			return i;
	}
	
	//we should never get here
	return -1;	
}

//look up face ID, used in print() for debugging
// FID = face ID
unsigned Mesh::FID(Face* f)
{
    for(unsigned i=0; i<faces.size(); i++)
	{
		if(faces[i] == f)
			return i;
	}
	
	//we should only get here if f is NULL (meaning the unbounded, outside face)
	return -1;	
}

//look up vertex ID, used in print() for debugging
// VID = vertex ID
unsigned Mesh::VID(Vertex* v)
{
    for(unsigned i=0; i<vertices.size(); i++)
    {
        if(vertices[i] == v)
            return i;
    }

    //we should only get here if f is NULL (meaning the unbounded, outside face)
    return -1;
}

//attempts to find inconsistencies in the DCEL arrangement
void Mesh::test_consistency()
{
    //check faces
    qDebug() << "Checking faces:";
    bool face_problem = false;
    std::set<int> edges_found_in_faces;

    for(std::vector<Face*>::iterator it = faces.begin(); it != faces.end(); ++it)
    {
        Face* face = *it;
        qDebug() << "  Checking face " << FID(face);

        if(face->get_boundary() == NULL)
        {
            qDebug() << "    PROBLEM: face" << FID(face) << "has null edge pointer.";
            face_problem = true;
        }
        else
        {
            Halfedge* start = face->get_boundary();
            edges_found_in_faces.insert(HID(start));

            if(start->get_face() != face)
            {
                qDebug() << "    PROBLEM: starting halfedge edge" << HID(start) << "of face" << FID(face) << "doesn't point back to face.";
                face_problem = true;
            }

            if(start->get_next() == NULL)
                qDebug() << "    PROBLEM: starting halfedge" << HID(start) << "of face" << FID(face) << "has NULL next pointer.";
            else
            {
                Halfedge* cur = start->get_next();
                int i = 0;
                while(cur != start)
                {
                    edges_found_in_faces.insert(HID(cur));

                    if(cur->get_face() != face)
                    {
                        qDebug() << "    PROBLEM: halfedge edge" << HID(cur) << "points to face" << FID(cur->get_face()) << "instead of face" << FID(face);
                        face_problem = true;
                        break;
                    }

                    if(cur->get_next() == NULL)
                    {
                        qDebug() << "    PROBLEM: halfedge" << HID(cur) << "has NULL next pointer.";
                        face_problem = true;
                        break;
                    }
                    else
                        cur = cur->get_next();

                    i++;
                    if(i >= 1000)
                    {
                        qDebug() << "    PROBLEM: halfedges of face" << FID(face) << "do not form a cycle (or, if they do, it has more than 1000 edges).";
                        face_problem = true;
                        break;
                    }
                }
            }

        }
    }//end face loop
    if(!face_problem)
        qDebug() << "   ---No problems detected among faces.";
    else
        qDebug() << "   ---Problems detected among faces.";

    //find exterior halfedges
    Halfedge* start = halfedges[1];
    Halfedge* cur = start;
    do{
        edges_found_in_faces.insert(HID(cur));

        if(cur->get_next() == NULL)
        {
            qDebug() << "    PROBLEM: halfedge " << HID(cur) << " has NULL next pointer.";
            break;
        }
        cur = cur->get_next();
    }while(cur != start);

    //check if all edges were found
    bool all_edges_found = true;
    for(unsigned i=0; i<halfedges.size(); i++)
    {
        if(edges_found_in_faces.find(i) == edges_found_in_faces.end())
        {
            qDebug() << "  PROBLEM: halfedge" << i << "not found in any face";
            all_edges_found = false;
        }
    }
    if(all_edges_found)
        qDebug() << "   ---All halfedges found in faces, as expected.";


    //check anchor lines
    qDebug() << "Checking anchor lines:\n";
    bool curve_problem = false;
    std::set<int> edges_found_in_curves;

/*    for(std::set<Anchor*>::iterator it = all_anchors.begin(); it != all_anchors.end(); ++it)
    {
        Anchor* anchor = *it;
        qDebug() << "  Checking line for anchor (" << anchor->get_x() <<"," << anchor->get_y() << ")";

        Halfedge* edge = anchor->get_line();
        do{
            edges_found_in_curves.insert(HID(edge));
            edges_found_in_curves.insert(HID(edge->get_twin()));

            if(edge->get_anchor() != anchor)
            {
                qDebug() << "    PROBLEM: halfedge" << HID(edge) << "does not point to this Anchor.";
                curve_problem = true;
            }
            if(edge->get_twin()->get_anchor() != anchor)
            {
                qDebug() << "    PROBLEM: halfedge" << HID(edge->get_twin()) << ", twin of halfedge " << HID(edge) << ", does not point to this anchor.";
                curve_problem = true;
            }

            if(edge->get_next() == NULL)
            {
                qDebug() << "    PROBLEM: halfedge" << HID(edge) << "has NULL next pointer.";
                curve_problem = true;
                break;
            }

            //find next edge in this line
            edge = edge->get_next();
            while(edge->get_anchor() != anchor)
                edge = edge->get_twin()->get_next();

        }while(edge->get_origin()->get_x() < INFTY);
    }//end anchor line loop
*/
    //ignore halfedges on both sides of boundary
    start = halfedges[1];
    cur = start;
    do{
        edges_found_in_curves.insert(HID(cur));
        edges_found_in_curves.insert(HID(cur->get_twin()));

        if(cur->get_next() == NULL)
        {
            qDebug() << "    PROBLEM: halfedge" << HID(cur) << "has NULL next pointer.";
            break;
        }
        cur = cur->get_next();
    }while(cur != start);

    //check if all edges were found
    all_edges_found = true;
    for(unsigned i=0; i<halfedges.size(); i++)
    {
        if(edges_found_in_curves.find(i) == edges_found_in_curves.end())
        {
            qDebug() << "  PROBLEM: halfedge" << i << "not found in any anchor line";
            all_edges_found = false;
        }
    }
    if(all_edges_found)
        qDebug() << "   ---All halfedges found in curves, as expected.";


    if(!curve_problem)
        qDebug() << "   ---No problems detected among anchor lines.";
    else
        qDebug() << "   ---Problems detected among anchor lines.";


    //check anchor lines
    qDebug() << "Checking order of vertices along right edge of the strip:";
    Halfedge* redge = halfedges[3];
    while(redge != halfedges[1])
    {
        qDebug() << " y = " << redge->get_origin()->get_y() << "at vertex" << VID(redge->get_origin());
        redge = redge->get_next();
    }
}//end test_consistency()



/********** the following objects and functions are for exact comparisons **********/

//Crossing constructor
//precondition: Anchors a and b must be comparable
Mesh::Crossing::Crossing(Anchor* a, Anchor* b, Mesh* m) : a(a), b(b), m(m)
{
    //store the x-coordinate of the crossing for fast (inexact) comparisons
    x = (m->y_grades[a->get_y()] - m->y_grades[b->get_y()])/(m->x_grades[a->get_x()] - m->x_grades[b->get_x()]);
}

//returns true iff this Crossing has (exactly) the same x-coordinate as other Crossing
bool Mesh::Crossing::x_equal(const Crossing *other) const
{
    if(Mesh::almost_equal(x, other->x)) //then compare exact values
    {
        //find exact x-values
        exact x1 = (m->y_exact[a->get_y()] - m->y_exact[b->get_y()])/(m->x_exact[a->get_x()] - m->x_exact[b->get_x()]);
        exact x2 = (m->y_exact[other->a->get_y()] - m->y_exact[other->b->get_y()])/(m->x_exact[other->a->get_x()] - m->x_exact[other->b->get_x()]);

        return x1 == x2;
    }
    //otherwise, x-coordinates are not equal
    return false;
}

//CrossingComparator for ordering crossings: first by x (left to right); for a given x, then by y (low to high)
bool Mesh::CrossingComparator::operator()(const Crossing* c1, const Crossing* c2) const	//returns true if c1 comes after c2
{
    //TESTING
    if(c1->a->get_position() >= c1->b->get_position() || c2->a->get_position() >= c2->b->get_position())
    {
        qDebug() << "INVERTED CROSSING ERROR\n";
        qDebug() << "crossing 1 involves anchors " << c1->a << " (pos " << c1->a->get_position() << ") and " << c1->b << " (pos " << c1->b->get_position() << "),";
        qDebug() << "crossing 2 involves anchors " << c2->a << " (pos " << c2->a->get_position() << ") and " << c2->b << " (pos " << c2->b->get_position() << "),";
        throw std::exception();
    }

    Mesh* m = c1->m;    //makes it easier to reference arrays in the mesh

    //now do the comparison
    //if the x-coordinates are nearly equal as double values, then compare exact values
    if(Mesh::almost_equal(c1->x, c2->x))
    {
        //find exact x-values
        exact x1 = (m->y_exact[c1->a->get_y()] - m->y_exact[c1->b->get_y()])/(m->x_exact[c1->a->get_x()] - m->x_exact[c1->b->get_x()]);
        exact x2 = (m->y_exact[c2->a->get_y()] - m->y_exact[c2->b->get_y()])/(m->x_exact[c2->a->get_x()] - m->x_exact[c2->b->get_x()]);

        //if the x-values are exactly equal, then consider the y-values
        if(x1 == x2)
        {
            //find the y-values
            double c1y = m->x_grades[c1->a->get_x()]*(c1->x) - m->y_grades[c1->a->get_y()];
            double c2y = m->x_grades[c2->a->get_x()]*(c2->x) - m->y_grades[c2->a->get_y()];

            //if the y-values are nearly equal as double values, then compare exact values
            if(Mesh::almost_equal(c1y, c2y))
            {
                //find exact y-values
                exact y1 = m->x_exact[c1->a->get_x()]*x1 - m->y_exact[c1->a->get_y()];
                exact y2 = m->x_exact[c2->a->get_x()]*x2 - m->y_exact[c2->a->get_y()];

                //if the y-values are exactly equal, then sort by relative position of the lines
                if(y1 == y2)
                    return c1->a->get_position() > c2->a->get_position();   //Is there a better way???

                //otherwise, the y-values are not equal
                return y1 > y2;
            }
            //otherwise, the y-values are not almost equal
            return c1y > c2y;
        }
        //otherwise, the x-values are not equal
        return x1 > x2;
    }
    //otherwise, the x-values are not almost equal
    return c1->x > c2->x;
}

//epsilon value for use in comparisons
double Mesh::epsilon = pow(2,-30);

//test whether two double values are almost equal (indicating that we should do exact comparison)
bool Mesh::almost_equal(const double a, const double b)
{
    double diff = std::abs(a - b);
    if(diff <= epsilon)
        return true;

    if(diff <= (std::abs(a) + std::abs(b))*epsilon)
        return true;
    return false;
}

