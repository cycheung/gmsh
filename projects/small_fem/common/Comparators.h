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

   @fn bool DofComparator::operator()(const Dof* a, const Dof* b) const;
   @param a A Dof pointer
   @param b A second Dof pointer (possibly pointing to the same Dof as @c a)
   @return Returns:
   @li @c true, if the Dof @em pointed by @c a is @em smaller 
   than the Dof @em pointed by @c b
   @li @c false, otherwise
 */


/**
   @class ElementComparator
   @brief A comparator for MElement%s

   This class is able to compare two MElement%s.

   @fn bool ElementComparator::operator()(const MElement* a, const MElement* b) const;
   @param a A MElement pointer
   @param b A second MElement pointer (possibly pointing to the same MElement as @c a)
   @return Returns:
   @li @c true, if the MElement @em pointed by @c a is @em smaller 
   than the MElement @em pointed by @c b
   @li @c false, otherwise
 */

/**
   @class VertexComparator
   @brief A comparator for MVertices

   This class is able to compare two MVertices.

   @fn bool VertexComparator::operator()(const MVertex* a, const MVertex* b) const;
   @param a A MVertex pointer
   @param b A second MVertex pointer (possibly pointing to the same MVertex as @c a)
   @return Returns:
   @li @c true, if the MVertex @em pointed by @c a is @em smaller 
   than the MVertex @em pointed by @c b
   @li @c false, otherwise
 */

/**
   @class EdgeComparator
   @brief A comparator for MEdge%s (without @em orientation notion)

   This class is able to compare two MEdge%s.

   @warning
   With this comparator, two MEdge%s with the @em same @em points,
   but @em different @em orientations, are handled as the @em same
   MEdge.

   @fn bool EdgeComparator::operator()(const MEdge* a, const MEdge* b) const;
   @param a A MEdge pointer
   @param b A second MEdge pointer (possibly pointing to the same MEdge as @c a)
   @return Returns:
   @li @c true, if the MEdge @em pointed by @c a is @em smaller 
   than the MEdge @em pointed by @c b
   @li @c false, otherwise

   @warning
   With this comparator, two MEdge%s with the @em same @em points,
   but @em different @em orientations, are handled as the @em same
   MEdge.
 */

/**
   @class FaceComparator
   @brief A comparator for MFace%s (without @em orientation notion)

   This class is able to compare two MFace%s.

   @warning
   With this comparator, two MFace%s with the @em same @em points,
   but @em different @em orientations, are handled as the @em same
   MFace.

   @fn bool FaceComparator::operator()(const MFace* a, const MFace* b) const;
   @param a A MFace pointer
   @param b A second MFace pointer (possibly pointing to the same MFace as @c a)
   @return Returns:
   @li @c true, if the MFace @em pointed by @c a is @em smaller 
   than the MFace @em pointed by @c b
   @li @c false, otherwise

   @warning
   With this comparator, two MFace%s with the @em same @em points,
   but @em different @em orientations, are handled as the @em same
   MFace.
 */

/**
   @class OrientedEdgeComparator
   @brief A comparator for MEdge%s (with @em orientation notion)

   This class is able to compare two MEdge%s.

   @warning
   With this comparator, two MEdge%s with the @em same @em points,
   but @em different @em orientations, are handled as @em different
   MEdge%s.

   @fn bool OrientedEdgeComparator::operator()(const MEdge* a, const MEdge* b) const;
   @param a A MEdge pointer
   @param b A second MEdge pointer (possibly pointing to the same MEdge as @c a)
   @return Returns:
   @li @c true, if the MEdge @em pointed by @c a is @em smaller 
   than the MEdge @em pointed by @c b
   @li @c false, otherwise

   @warning
   With this comparator, two MEdge%s with the @em same @em points,
   but @em different @em orientations, are handled as @em different
   MEdge%s.
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
