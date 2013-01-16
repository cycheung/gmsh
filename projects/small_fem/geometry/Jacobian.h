#ifndef _JACOBIAN_H_
#define _JACOBIAN_H_

#include <vector>
#include <map>

#include "Comparators.h"
#include "GroupOfElement.h"
#include "fullMatrix.h"

/**
   @class Jacobian
   @brief Handels Jacobians of a Group of Element

   This class handels the Jacobians of a
   Group of Element (GroupOfElement).@n

   The Jacobian%s will be computed at a given point.@n
   This point is given at construction time.@n

   For this class the jacobian matrix is defined as:@n@n
   @f$J~=~\left(
   \begin{array}{ccc}
   \displaystyle\frac{\partial{}x}{\partial{}u} &
   \displaystyle\frac{\partial{}y}{\partial{}u} &
   \displaystyle\frac{\partial{}z}{\partial{}u}\\
   \displaystyle\frac{\partial{}x}{\partial{}v} &
   \displaystyle\frac{\partial{}y}{\partial{}v} &
   \displaystyle\frac{\partial{}z}{\partial{}v}\\
   \displaystyle\frac{\partial{}x}{\partial{}w} &
   \displaystyle\frac{\partial{}y}{\partial{}w} &
   \displaystyle\frac{\partial{}z}{\partial{}w}\\
   \end{array}
   \right)@f$
*/

class Jacobian{
 private:
  const std::vector<const MElement*>* element;
  const fullMatrix<double>*           point;
  unsigned int                        row;

  typedef std::pair<const fullMatrix<double>*, double> jac_t;

  std::map<const MElement*, jac_t*, ElementComparator>* jac;
  std::map<const MElement*, jac_t*, ElementComparator>* invJac;

 public:
  Jacobian(const GroupOfElement& goe,
	   const fullMatrix<double>& point,
	   unsigned int row);

  ~Jacobian(void);

  void computeJacobians(void);
  void computeInvertJacobians(void);

  const std::pair<const fullMatrix<double>*, double>&
    getJacobian(const MElement& element) const;

  const std::pair<const fullMatrix<double>*, double>&
    getInvertJacobian(const MElement& element) const;

  const std::vector<const MElement*>&
    getAllElements(void) const;

  fullVector<double> getPoint(void) const;

 private:
  void deleteJac(void);
  void deleteInvertJac(void);

  void computeInvertFromJac(void);
  void computeInvertFromScratch(void);
};

/**
   @fn Jacobian::Jacobian
   @param goe A Group of Element (GroupOfElement)
   @param point A @c [ @c N @c x @c 3 @c ] matrix
   (a set of @c N points coordinates)
   @param row A row of the matrix @c point
   (row is ranging from @c 0 to @c N @c - @c 1)

   Instanciates a new Jacobian@n

   The given group of element defines the elements
   where jacobians will be computed@n
   @see Jacobian::getAllElements()

   The given matrix and its row defines the point
   where jacobians will be computed@n
   @see Jacobian::getPoint()
   **

   @fn Jacobian::~Jacobian

   Deletes this Jacobian
   **

   @fn Jacobian::computeJacobians

   Computes the jacobians (and its determinant)
   of Jacobian::getAllElements() at Jacobian::getPoint()
   **

   @fn Jacobian::computeInvertJacobians

   Computes the @em invert of jacobians
   (and the @em non @em invert determinant)
   of Jacobian::getAllElements() at Jacobian::getPoint()
   **

   @fn Jacobian::getJacobian
   @param element A MElement
   @return Returns a pair with the jacobian matrix,
   and its determinant, of Jacobian::getAllElements()
   evaluated at Jacobian::getPoint()

   @note
   If no jacobian has been precomputed,
   an Exception is thrown

   @see Jacobian::computeJacobians()
   **

   @fn Jacobian::getInvertJacobian
   @param element A MElement
   @return Returns a pair with the @em invert
   jacobian matrix, and its @em non @em invert
   determinant, of Jacobian::getAllElements()
   evaluated at Jacobian::getPoint()

   @note
   If no jacobian has been precomputed,
   an Exception is thrown

   @see Jacobian::computeJacobians()
   **

   @fn Jacobian::getAllElements
   @return Returns the elements where the jacobians
   will be computed
   **

   @fn Jacobian::getPoint
   @return Returns the point where the jacobians
   will be computed
 */


//////////////////////
// Inline Functions //
//////////////////////

inline const std::vector<const MElement*>&
Jacobian::getAllElements(void) const{
  return *element;
}

inline fullVector<double>Jacobian::
getPoint(void) const{
  fullVector<double> p(3);

  p(0) = (*point)(row, 0);
  p(1) = (*point)(row, 1);
  p(2) = (*point)(row, 2);

  return p;
}

inline void Jacobian::
computeInvertJacobians(void){
  if(invJac)
    return;

  else if(jac)
    computeInvertFromJac();

  else
    computeInvertFromScratch();
}

#endif
