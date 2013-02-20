#ifndef _ORIENTATIONSORT_H_
#define _ORIENTATIONSORT_H_

#include <vector>

#include "MElement.h"
#include "Basis.h"

/**
   @class OrientationSort
   @brief Sorting Predicate for Elements orientation

   Sorting Predicate for Elements
   with respect to Basis::getOrientation()
 */

class OrientationSort{
 private:
  const Basis* basis;

 public:
   OrientationSort(const Basis& basis);
  ~OrientationSort(void);

  const Basis& getBasis(void) const;

  bool operator()(const MElement* a, const MElement* b) const;
};

/**
   @fn OrientationSort::OrientationSort
   @param basis A Basis
   Instanciates a new OrientationSort,
   with the given Basis for sorting predicate
   @see Basis::getOrientation()
   **

   @fn OrientationSort::~OrientationSort
   Deletes this OrientationSort
   **

   @fn OrientationSort::getBasis
   @return Returns the Basis used for orienting
   elements
   @see OrientationSort::OrientationSort
   **

   @fn OrientationSort::operator()(const MElement*, const MElement*) const
   @param a A MElement
   @param b A MElement (possibly the same as @c a)
   @return Returns:
   @li @c true, if
   getBasis().getOrientation(*a) < getBasis().getOrientation(*b)
   @li @c false, otherwise
 */


//////////////////////
// Inline Functions //
//////////////////////

inline const Basis& OrientationSort::getBasis(void) const{
  return *basis;
}

inline bool OrientationSort::operator()
(const MElement* a, const MElement* b) const{
  return
    basis->getOrientation(*a) < basis->getOrientation(*b);
}

#endif
