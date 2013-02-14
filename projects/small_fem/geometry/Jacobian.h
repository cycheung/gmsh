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

   The Jacobian%s will be computed at a given set of points.@n

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

   @todo
   GMSH WARNING: Memory leak with MElement::getJacobian() !!
*/

class Jacobian{
 private:
  const GroupOfElement*     goe;
  const fullMatrix<double>* point;
  unsigned int              nElement;
  unsigned int              nPoint;

  const std::vector<std::pair<const MElement*, ElementData> >* element;

  typedef std::pair<const fullMatrix<double>*, double> jac_pair;
  typedef std::vector<const jac_pair*>                 jac_t;

  std::map<const MElement*, jac_t*, ElementComparator>* jac;
  std::map<const MElement*, jac_t*, ElementComparator>* invJac;

 public:
  Jacobian(const GroupOfElement& goe,
	   const fullMatrix<double>& point);

  ~Jacobian(void);

  void computeJacobians(void);
  void computeInvertJacobians(void);

  const std::vector<const std::pair<const fullMatrix<double>*, double>*>&
    getJacobian(const MElement& element) const;

  const std::vector<const std::pair<const fullMatrix<double>*, double>*>&
    getInvertJacobian(const MElement& element) const;

  const GroupOfElement&     getAllElements(void) const;
  const fullMatrix<double>& getAllPoints(void) const;

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
   (a set of @c N points and their coordinates)

   Instanciates a new Jacobian@n

   The given group of element defines the elements
   where jacobians will be computed@n
   @see Jacobian::getAllElements()

   The given matrix defines the points
   where jacobians will be computed@n
   @see Jacobian::getAllPoints()
   **

   @fn Jacobian::~Jacobian

   Deletes this Jacobian
   **

   @fn Jacobian::computeJacobians

   Computes the jacobians (and its determinant)
   of all Jacobian::getAllElements()
   at all Jacobian::getAllPoints()
   **

   @fn Jacobian::computeInvertJacobians

   Computes the @em invert of jacobians
   (and the @em non @em invert determinant)
   of all Jacobian::getAllElements()
   at all Jacobian::getAllPoints()
   **

   @fn Jacobian::getJacobian
   @param element A MElement
   @return Returns a vector of pairs.@n
   The @c i-th element of this vector is such that:

   @li Its first entry is the jacobian matrix
   of the given element
   evaluated at Jacobian::getAllPoints()(i, :)

   @li Its second entry is the jacobian matrix
   detetminant of the given element
   evaluated at Jacobian::getAllPoints()(i, :)

   @note
   If no jacobian has been precomputed,
   an Exception is thrown

   @see Jacobian::computeJacobians()
   **

   @fn Jacobian::getInvertJacobian
   @param element A MElement
   @return Returns a vector of pairs.@n
   The @c i-th element of this vector is such that:

   @li Its first entry is the @em invert
   jacobian matrix of the given element
   evaluated at Jacobian::getAllPoints()(i, :)

   @li Its second entry is the
   @em non @em inverted jacobian matrix
   detetminant of the given element
   evaluated at Jacobian::getAllPoints()(i, :)

   @note
   If no inverted jacobian has been precomputed,
   an Exception is thrown

   @see Jacobian::computeInvertJacobians()
   **

   @fn Jacobian::getAllElements
   @return Returns the Group Of Elements
   (see GroupOfElement) where the jacobians
   will be computed
   **

   @fn Jacobian::getAllPoints
   @return Returns the points where the jacobians
   will be computed
 */


//////////////////////
// Inline Functions //
//////////////////////

inline const GroupOfElement&
Jacobian::getAllElements(void) const{
  return *goe;
}

inline const fullMatrix<double>& Jacobian::
getAllPoints(void) const{
  return *point;
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
