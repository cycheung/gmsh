#ifndef _GROUPOFJACOBIAN_H_
#define _GROUPOFJACOBIAN_H_

#include "GroupOfElement.h"
#include "fullMatrix.h"
#include "Jacobian.h"

#include <string>

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
                  const fullMatrix<double>& point,
                  const std::string type);

  ~GroupOfJacobian(void);

  const Jacobian&       getJacobian(unsigned int i) const;
  const GroupOfElement& getAllElements(void) const;
};

/////////////////////
// Inline Function //
/////////////////////

inline const Jacobian&
GroupOfJacobian::getJacobian(unsigned int i) const{
  return *jac[i];
}

inline const GroupOfElement&
GroupOfJacobian::getAllElements(void) const{
  return *goe;
}

#endif
