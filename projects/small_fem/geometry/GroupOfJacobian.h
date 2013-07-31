#ifndef _GROUPOFJACOBIAN_H_
#define _GROUPOFJACOBIAN_H_

#include <string>

#include "GroupOfElement.h"
#include "fullMatrix.h"
#include "Basis.h"

#include "Jacobian.h"

/**
   @class GroupOfJacobian
   @brief A set of Jacobian%s

   A set of Jacobian%s defined over a GroupOfElement.
 */

class GroupOfJacobian{
 private:
  const GroupOfElement* goe;
  Jacobian**            jac;

 public:
  GroupOfJacobian(const GroupOfElement& goe,
                  const Basis& basis,
                  const fullMatrix<double>& point,
                  const std::string type);

  ~GroupOfJacobian(void);

  const Jacobian&       getJacobian(size_t i) const;
  const GroupOfElement& getAllElements(void) const;

  std::string toString(void) const;
};

/**
   @fn GroupOfJacobian::GroupOfJacobian
   @param goe A group of element (GroupOfElement)
   @param basis A Basis
   @param point A [ N x 3 ] matrix (a set of N points and their coordinates)
   @param type A string

   Instanciate a GroupOfJacobian. This group will contain the Jacobian
   of the element of the given GroupOfElement, computed at the given points,
   for the given type, and with given Basis ReferenceSpace.

   @see Jacobian::Jacobian(const MElement&, const Basis&, const fullMatrix<double>&, const std::string)
   **

   @fn GroupOfJacobian::~GroupOfJacobian
   Deletes this GroupOfJacobian
   **

   @fn GroupOfJacobian::getJacobian
   @param i An integer
   @return Returns the Jacobian assiciated to the i-th
   Element of this GroupOfJacobian
   (and the i-th Element of the GroupOfElement given at instanciation)

   @see GroupOfJacobian::getAllElements()
   @see GroupOfElement::get();
   **

   @fn GroupOfJacobian::getAllElements
   @return Returns the GroupOfElement given at instanciation
   **

   @fn GroupOfJacobian::toString
   @return Returns a string describing this GroupOfJacobian
   **
 */

/////////////////////
// Inline Function //
/////////////////////

inline const Jacobian&
GroupOfJacobian::getJacobian(size_t i) const{
  return *jac[i];
}

inline const GroupOfElement&
GroupOfJacobian::getAllElements(void) const{
  return *goe;
}

#endif
