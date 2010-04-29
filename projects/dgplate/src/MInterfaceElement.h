//
// C++ Interface: terms
//
// Description: Class of interface element used for DG 2D only for the moment
// thus the interface element is a line
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
# ifndef _MINTERFACEELEMENT_H_
# define _MINTERFACEELEMENT_H_

#include "MLine.h"
#include "SVector3.h"
#include "MEdge.h"
#include "quadratureRules.h"
// closure + ghost cell (//) c'est une BC
class MInterfaceElement : public MLineN{ // or don't derivate but in this case a vector with the vertices of interface element has to be save ??
  protected :
    // table of pointer on the two elements linked to the interface element
    MElement *_numElem[2];
    // edge's number linked to interface element of minus and plus element
    int _numEdge[2];
    // dir = true if the edge and the interface element are defined in the same sens and dir = false otherwise
    bool _dir[2];
    // compute the characteritic size of one element (This function can be defined as a method of MElement) ??
    double characSize(MElement *e);

  public :
    MInterfaceElement(std::vector<MVertex*> &v, int num = 0, int part = 0, MElement *e_minus = 0, MElement *e_plus = 0);

    // Destructor
    ~MInterfaceElement(){}

    // Give the number of the elements in a vector
    MElement ** getElem(){return _numElem;}

    // Give the number of minus 0 or plus 1 element
    MElement * getElem(int index){return _numElem[index];}
    // Number of InterfaceElement == number of minus element if no int is given (cheat for dof on interface)
    int getNum() const{return _numElem[0]->getNum();};

    void getuvOnElem(const double u, double &uem, double &vem, double &uep, double &vep);

    // Return the edge number of element
    int getEdgeNumber(const int i) const {return _numEdge[i];}

    // Return the local vertex number of interface
    void getLocalVertexNum(const int i,std::vector<int> &vn);
    // Compute the characteristic size of the side h_s = max_e (area_e/perimeter_e)
    double getCharacteristicSize(){
      double cm = this->characSize(_numElem[0]);
      double cp = this->characSize(_numElem[1]);
      return cm > cp ? cm : cp;
    }
    bool getdir(const int i){return _dir[i];}
};

// Class used to build 2D interface element
class Iedge{
  protected :
    std::vector<MVertex*> vertex;
    MElement *elem;
    int phys;

  public :
    Iedge(std::vector<MVertex*> v,MElement *e, int i) : vertex(v), elem(e), phys(i){};
    ~Iedge(){};
    unsigned long int getkey(){
      int i1,i2,i3;
      i1 = vertex[0]->getNum();
      i2 = vertex[1]->getNum();
      i1>i2 ? i3=i1*100000+i2 : i3=i2*100000+i1; // change this
      return i3;
    }
    std::vector<MVertex*> getVertices() const{return vertex;}
    MElement* getElement() const{return elem;}
    MVertex* getFirstInteriorVertex() const {return vertex[2];}
    int getPhys() const {return phys;}
};
#endif // MInterfaceElement
