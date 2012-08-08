#ifndef _GROUPOFEDGE_H_
#define _GROUPOFEDGE_H_

#include <string>
#include <vector>
#include <map>

#include "Mesh.h"
#include "Group.h"
#include "GroupOfElement.h"
#include "MEdge.h"

/**
   @class GroupOfEdge
   @brief A Group of MEdge

   This class is collection of @em Edges (MEdge).@n
   This class is @em Group.
*/

class Mesh;
class GroupOfElement;

class GroupOfEdge: public Group{
 private:
  class EdgeComparator{
  public:
    bool operator()(const MEdge& e1, const MEdge& e2) const;
  };

 private:
  Mesh* mesh;

  static unsigned int nextId;
  unsigned int            id;

  unsigned int         nEdge;
  
  std::vector<MEdge*>*  edge;
  std::map<MEdge, int, EdgeComparator>* orientation;

 public:
  GroupOfEdge(const GroupOfElement& goe, Mesh& mesh);
  virtual ~GroupOfEdge(void);

  virtual unsigned int getNumber(void) const;
  virtual unsigned int getId(void) const;
  virtual unsigned int getType(void)   const;

  MEdge&                     get(unsigned int i) const;  
  const std::vector<MEdge*>& getAll(void) const;  
  int                        getOrientation(const MEdge& edge) const;
  Mesh&                      getMesh(void) const;

  virtual std::string toString(void) const;

 private:
  static MEdge invert(MEdge& edge);
};


//////////////////////
// Inline Functions //
//////////////////////

inline unsigned int GroupOfEdge::getNumber(void) const{
  return nEdge;
}

inline unsigned int GroupOfEdge::getId(void) const{
  return id;
}

inline unsigned int GroupOfEdge::getType(void) const{
  return 3;
}

inline MEdge& GroupOfEdge::get(unsigned int i) const{
  return (*(*edge)[i]);
}

inline const std::vector<MEdge*>& 
GroupOfEdge::getAll(void) const{
  return *edge;
}

inline Mesh& GroupOfEdge::getMesh(void) const{
  return *mesh;
}

inline bool GroupOfEdge::EdgeComparator::
operator()(const MEdge& e1, const MEdge& e2) const{
  return 
    ( e1.getVertex(0)->getNum() <  e2.getVertex(0)->getNum()) ||
    ((e1.getVertex(0)->getNum() == e2.getVertex(0)->getNum()) && 
     (e1.getVertex(1)->getNum() <  e2.getVertex(1)->getNum()));
}


#endif
