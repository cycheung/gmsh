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

   A set of Jacobian%s
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
