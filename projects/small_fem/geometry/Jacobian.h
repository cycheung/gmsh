#ifndef _JACOBIAN_H_
#define _JACOBIAN_H_

#include <vector>
#include <string>

#include "MElement.h"
#include "Exception.h"
#include "fullMatrix.h"

/**
   @class Jacobian
   @brief Handels Jacobians of an Element

   This class handels the Jacobians
   of an Element (MElement).@n

   The Jacobian%s will be computed at a given set
   of points.@n

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
  static const std::string jacobianString;
  static const std::string invertString;
  static const std::string bothString;

 private:
  typedef std::pair<const fullMatrix<double>*, double> jac_pair;
  typedef std::vector<const jac_pair*>                 jac_t;

 private:
  const MElement*         element;
  const fullMatrix<double>* point;
  std::string                type;

  jac_t* jac;
  jac_t* invJac;

 public:
  Jacobian(const MElement& element,
	   const fullMatrix<double>& point,
           const std::string type);

  ~Jacobian(void);

  const std::vector<const std::pair<const fullMatrix<double>*, double>*>&
    getJacobianMatrix(void) const;

  const std::vector<const std::pair<const fullMatrix<double>*, double>*>&
    getInvertJacobianMatrix(void) const;

  const MElement&           getElement(void) const;
  const fullMatrix<double>& getPoints(void) const;
  const std::string&        getType(void) const;

 private:
  void deleteJac(void);
  void deleteInvertJac(void);

  void computeJacobians(void);
  void computeInvertFromJac(void);
  void computeInvertFromScratch(void);
};

/**
   @fn Jacobian::Jacobian
   @param element An Element (MElement)
   @param point A @c [ @c N @c x @c 3 @c ] matrix
   (a set of @c N points and their coordinates)
   @param type A string

   Instanciates a new Jacobian of the given @em type:@n

   @li @c jacobian, to compute the jacobian matrices
   @li @c invert, to comupte the @em inverted jacobian matrices
   @li @c both, to compute @em both inverted and non inverted
   jacobian matrices@n@n

   The given matrix defines the set of points
   where jacobian matrices will be computed@n
   **

   @fn Jacobian::~Jacobian

   Deletes this Jacobian
   **

   @fn Jacobian::getJacobianMatrix
   @return Returns a vector of pairs.@n
   The @c i-th element of this vector is such that:

   @li Its first entry is the jacobian matrix
   evaluated at Jacobian::getPoints(i, :)

   @li Its second entry is the jacobian matrix
   determinant evaluated at Jacobian::getPoints(i, :)

   @note
   If Jacobian::getType() is @c invert,
   this method throws an Exception
   **

   @fn Jacobian::getInvertJacobian
   @return Returns a vector of pairs.@n
   The @c i-th element of this vector is such that:

   @li Its first entry is the @em invert
   jacobian matrix evaluated at Jacobian::getPoints(i, :)

   @li Its second entry is the
   @em non @em inverted jacobian matrix
   detetminant evaluated at Jacobian::getPoints(i, :)

   @note
   If the Jacobian::getType() is @c jacobian,
   this method throws an Exception
   **

   @fn Jacobian::getElement
   @return Return the Element on which
   the jacobians were computed
   **

   @fn Jacobian::getPoints
   @return Return the set of points on which
   the jacobians were computed
   **

   @fn Jacobian::getType
   @return Return the type of jacobians
   that were computed:

   @li @c jacobian, for the jacobian matrices
   @li @c invert, for the @em inverted jacobian matrices
   @li @c both, for @em both inverted and non inverted
   jacobian matrices
 */


//////////////////////
// Inline Functions //
//////////////////////

inline const std::vector<const std::pair<const fullMatrix<double>*, double>*>&
  Jacobian::getJacobianMatrix(void) const{

  if(jac)
    return *jac;

  else
    throw Exception
      ("This Jacobian type is %s -- can't get non inverted jacobian",
       type.c_str());
}

inline const std::vector<const std::pair<const fullMatrix<double>*, double>*>&
  Jacobian::getInvertJacobianMatrix(void) const{

  if(invJac)
    return *invJac;

  else
    throw Exception
      ("This Jacobian type is %s -- can't get inverted jacobian",
       type.c_str());
}

inline const MElement& Jacobian::getElement(void) const{
  return *element;
}

inline const fullMatrix<double>& Jacobian::getPoints(void) const{
  return *point;
}

inline const std::string& Jacobian::getType(void) const{
  return type;
}

#endif
