
#ifndef __ANCHOR_H__
#define __ANCHOR_H__

//forward declarations
class Halfedge;
struct xiMatrixEntry;


/**
 * \class	Anchor
 * \brief	Stores an Anchor: a multi-index pair along with a pointer to the line representing the Anchor in the arrangement
 * \author	Matthew L. Wright
 * \date	March 2014
 */
class Anchor
{
    public:
        Anchor(unsigned x, unsigned y);     //constructor, requires only x- and y-coordinates
        Anchor(const Anchor& other);        //copy constructor

        Anchor& operator= (const Anchor& other);	//assignment operator
        bool operator== (const Anchor& other) const;	//equality operator

        bool comparable(Anchor* other) const;   //tests whether two Anchors are (strongly) comparable

        unsigned get_x() const;		//get the discrete x-coordinate
        unsigned get_y() const;		//get the discrete y-coordinate

        void set_line(Halfedge* e);	//set the pointer to the line corresponding to this Anchor in the arrangement
        Halfedge* get_line() const;		//get the pointer to the line corresponding to this Anchor in the arrangement

        void set_position(unsigned p);  //sets the relative position of the Anchor line at the sweep line, used for Bentley-Ottmann DCEL construction algorithm
        unsigned get_position() const;  //gets the relative position of the Anchor line at the sweep line, used for Bentley-Ottmann DCEL construction algorithm

        bool is_above();       //returns true iff this Anchor is above the current slice line, used for the vineyard-update process of storing persistence data in cells of the arrangement
        void toggle();         //toggles above/below state of this Anchor; called whever the slice line crosses this Anchor in the vineyard-update process of storing persistence data

        xiMatrixEntry* get_entry(); //accessor

        void set_weight(unsigned long w);   //sets the estimate of the cost of updating the RU-decomposition when crossing this anchor
        unsigned long get_weight();         //returns estimate of the cost of updating the RU-decomposition when crossing this anchor

    private:
        unsigned x_coord;	//discrete x-coordinate
        unsigned y_coord;	//discrete y-coordinate

        xiMatrixEntry* entry;   //xiMatrixEntry at the position of this anchor

        Halfedge* dual_line;    //pointer to left-most halfedge corresponding to this Anchor in the arrangement
        unsigned position;      //relative position of Anchor line at sweep line, used for Bentley-Ottmann DCEL construction algorithm
        bool above_line;        //true iff this Anchor is above the current slice line, used for the vineyard-update process of storing persistence data in cells of the arrangement
        unsigned long weight;   //estimate of the cost of updating the RU-decomposition when crossing this anchor
};



class Anchor_LeftComparator_Full    //compares anchors to determine ordering at the x=-infty position
{
    public:
        bool operator() (const Anchor* lhs, const Anchor* rhs) const  //returns true if lhs comes before rhs
        {
            if(lhs->get_x() > rhs->get_x())		//first compare x-coordinates (descending)
                return true;
            if(lhs->get_x() == rhs->get_x() && lhs->get_y() > rhs->get_y())		//then compare y-coordinates (descending)
                return true;
            return false;
        }
};


class Anchor_LeftComparator_Half    //compares anchors to determine ordering at the x=0 position
{
    public:
        bool operator() (const Anchor* lhs, const Anchor* rhs) const  //returns true if lhs comes before rhs
        {
            if(lhs->get_y() > rhs->get_y())		//first compare y-coordinates (reverse order)
                return true;
            if(lhs->get_y() == rhs->get_y() && lhs->get_x() < rhs->get_x())		//then compare x-coordinates (natural order)
                return true;
            return false;
        }
};

#endif // __ANCHOR_H__
