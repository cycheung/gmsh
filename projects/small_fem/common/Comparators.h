#ifndef _COMPARATORS_H_
#define _COMPARATORS_H_

#include "Dof.h"
#include "MElement.h"
#include "MVertex.h"
#include "MEdge.h"

/**
   @class DofComparator
   @brief A comparator for Dof%s

   This class is able to compare two Dof%s.

   @fn bool DofComparator::operator()(const Dof*, const Dof*)
   @param a A Dof pointer
   @param b A second Dof pointer
   @return Returns:
   @li true, if the Dof pointed by a is smaller
   than the Dof pointed by b
   @li false, otherwise
 */


/**
   @class ElementComparator
   @brief A comparator for MElement%s

   This class is able to compare two MElement%s.

   @fn bool ElementComparator::operator()(const MElement*, const MElement*)
   @param a A MElement pointer
   @param b A second MElement pointer
   @return Returns:
   @li true, if the MElement pointed by a is smaller
   than the MElement pointed by b
   @li false, otherwise
 */

/**
   @class VertexComparator
   @brief A comparator for MVertices

   This class is able to compare two MVertices.

   @fn bool VertexComparator::operator()(const MVertex*, const MVertex*)
   @param a A MVertex pointer
   @param b A second MVertex pointer
   @return Returns:
   @li true, if the MVertex pointed by a is smaller
   than the MVertex pointed by b
   @li false, otherwise
 */

/**
   @class EdgeComparator
   @brief A comparator for MEdge%s (without orientation notion)

   This class is able to compare two MEdge%s (without orientation notion).

   With this comparator, two MEdge%s with the same points,
   but different orientations, are handled as the same MEdge.

   @fn bool EdgeComparator::operator()(const MEdge*, const MEdge*) const
   @param a A MEdge pointer
   @param b A second MEdge pointer
   @return Returns:
   @li true, if the MEdge pointed by a is smaller
   than the MEdge pointed by b
   @li false, otherwise
 */

/**
   @class FaceComparator
   @brief A comparator for MFace%s (without orientation notion)

   This class is able to compare two MFace%s (without orientation notion).

   With this comparator, two MFace%s with the same points,
   but different orientations, are handled as the same MFace.

   @fn bool FaceComparator::operator()(const MFace*, const MFace*)
   @param a A MFace pointer
   @param b A second MFace pointer
   @return Returns:
   @li true, if the MFace pointed by a is smaller
   than the MFace pointed by b
   @li false, otherwise
 */

/**
   @class OrientedEdgeComparator
   @brief A comparator for MEdge%s (with orientation notion)

   This class is able to compare two MEdge%s (with orientation notion).

   With this comparator, two MEdge%s with the same points,
   but different orientations, are handled as different MEdge%s.

   @fn bool OrientedEdgeComparator::operator()(const MEdge*, const MEdge*)
   @param a A MEdge pointer
   @param b A second MEdge pointer
   @return Returns:
   @li true, if the MEdge pointed by a is smaller
   than the MEdge pointed by b
   @li false, otherwise
 */

class DofComparator{
 public:
  bool operator()(const Dof* a, const Dof* b) const;
};

class ElementComparator{
 public:
  bool operator()(const MElement *a, const MElement *b) const;
};

class VertexComparator{
 public:
  bool operator()(const MVertex *a, const MVertex *b) const;
};

class EdgeComparator{
 public:
  bool operator()(const MEdge *a, const MEdge *b) const;
};

class FaceComparator{
 public:
  bool operator()(const MFace *a, const MFace *b) const;
};

class OrientedEdgeComparator{
 public:
  bool operator()(const MEdge* a, const MEdge* b) const;
};


//////////////////////
// Inline Functions //
//////////////////////

inline bool DofComparator::
operator()(const Dof* a, const Dof* b) const{
  return *a < *b;
}

inline bool ElementComparator::
operator()(const MElement *a, const MElement *b) const{
  return a->getNum() < b->getNum();
}

inline bool VertexComparator::
operator()(const MVertex *a, const MVertex *b) const{
  return a->getNum() < b->getNum();
}

inline bool EdgeComparator::
operator()(const MEdge *a, const MEdge *b) const{
  if(a->getMinVertex() < b->getMinVertex()) return true;
  if(a->getMinVertex() > b->getMinVertex()) return false;
  if(a->getMaxVertex() < b->getMaxVertex()) return true;
  return false;
}

inline bool OrientedEdgeComparator::
operator()(const MEdge* a, const MEdge* b) const{
  return
    ( a->getVertex(0)->getNum() <  b->getVertex(0)->getNum()) ||
    ((a->getVertex(0)->getNum() == b->getVertex(0)->getNum()) &&
     (a->getVertex(1)->getNum() <  b->getVertex(1)->getNum()));
}


#endif
