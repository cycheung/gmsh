#ifndef _COMPARATORS_H_
#define _COMPARATORS_H_

#include "Group.h"
#include "MElement.h"
#include "MEdge.h"

class GroupComparator{
 public:
  bool operator()(const Group* a, const Group* b) const;
};

class ElementComparator{
 public:
  bool operator()(const MElement *e1, const MElement *e2) const;
};

class EdgeComparator{
 public:
  bool operator()(const MEdge *e1, const MEdge *e2) const;
};

class OrientedEdgeComparator{
 public:
  bool operator()(const MEdge* e1, const MEdge* e2) const;
};


//////////////////////
// Inline Functions //
//////////////////////

inline bool GroupComparator::
operator()(const Group* a, const Group* b) const{
  return a->getId() < b->getId();
}

inline bool ElementComparator::
operator()(const MElement *e1, const MElement *e2) const{
  return e1->getNum() < e2->getNum();
}

inline bool EdgeComparator::
operator()(const MEdge *e1, const MEdge *e2) const{
  if(e1->getMinVertex() < e2->getMinVertex()) return true;
  if(e1->getMinVertex() > e2->getMinVertex()) return false;
  if(e1->getMaxVertex() < e2->getMaxVertex()) return true;
  return false;
}

inline bool OrientedEdgeComparator::
operator()(const MEdge* e1, const MEdge* e2) const{
  return 
    ( e1->getVertex(0)->getNum() <  e2->getVertex(0)->getNum()) ||
    ((e1->getVertex(0)->getNum() == e2->getVertex(0)->getNum()) && 
     (e1->getVertex(1)->getNum() <  e2->getVertex(1)->getNum()));
}


#endif
