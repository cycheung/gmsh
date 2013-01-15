#ifndef _JACOBIAN_H_
#define _JACOBIAN_H_

#include <map>

#include "Comparators.h"
#include "GroupOfElement.h"
#include "fullMatrix.h"

/**
   @class Jacobian
   @brief Handels Jacobians of a Group of Element

   This class handels the Jacobians of a 
   Group of Element (GroupOfElement)
*/

class Jacobian{
 private:
  typedef std::pair<fullMatrix<double>*, double> jac_t;

  std::map<const MElement*, jac_t*, ElementComparator>* jac;

 public:
  Jacobian(const GroupOfElement& goe, 
	   const fullMatrix<double>& point,
	   unsigned int row);
  
  ~Jacobian(void);
  
  
  
};

/**
   @fn Jacobian::Jacobian
   @param goe A Group of Element (GroupOfElement)
   @param point A @c [ @c N @c x @c 3 @c ] matrix
   (a set of @c N points coordinates)
   @param row A row of the matrix @c point
   (row is ranging from @c 0 to @c N @c - @c 1)

   Instanciates a new Jacobian@n
   
   The @em jacobian @em matrix 
   (and the determinant) of every Element
   of the given GroupOfElement will be computed
   at the given point (the given row of @c point),
   in the Element @em Reference @em Space
   **

   @fn Jacobian::~Jacobian
   
   Deletes this Jacobian
 */

#endif
