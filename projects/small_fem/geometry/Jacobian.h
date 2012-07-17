#ifndef _JACOBIAN_H_
#define _JACOBIAN_H_

#include <vector>
#include "fullMatrix.h"
#include "MElement.h"

/**
   @class Jacobian
   @brief Transformations between
   physical and reference spaces

   This class handles the transformations
   between physical and reference spaces.@n

   The @em pysical space is defined by:
   @li @c X and @c Y coordinates in @c 2D
   @li @c X, @c Y and @c Z coordinates in @c 3D

   The @em reference space is defined by:
   @li @c U and @c V coordinates in @c 2D
   @li @c U, @c V and @c W coordinates in @c 3D

   @todo
   Implement a cache ?@n
   const_cast is dirty:
   Change MElement with const qualifier
*/

class Jacobian{
 private:
  MElement* element;

 public:
   Jacobian(const MElement& element);
  ~Jacobian(void);

  double det(void) const;
  double det(const fullVector<double>& UVW) const;

  fullVector<double> map(const fullVector<double>& UVW) const;

  fullVector<double> grad(const fullVector<double>& gradUVW) const;

  fullVector<double> invMap(const fullVector<double>& XYZ) const;
};

/**
   @fn Jacobian::Jacobian
   @param element The Element from which the Jacobian 
   will be computed
   @return Returns a new Jacobian

   @fn Jacobian::~Jacobian
   @return Deletes this Jacobian

   @fn Jacobian::det(const fullVector<double>&) const
   @param UVW A @c 3D Vector with the coordinate
   of a point in the @em reference (@c 3D) space
   @return Returns the determinant of the 
   jacobian matrix at the given point

   @fn Jacobian::map(const fullVector<double>&) const
   @param UVW A @c 3D Vector with the coordinate 
   of a point in the @em reference (@c 3D) space
   @returns Returns the coordiantes of the given point
   in the @em physical space

   @fn Jacobian::grad(const fullVector<double>&) const
   @param gradUVW A gradient in the @em reference space
   @returns Returns the given gradient in the 
   @em physical space

   @fn Jacobian::invMap(const fullVector<double>&) const
   @param XYZ A @c 3D Vector with the coordinate 
   of a point in the @em physical (@c 3D) space
   @returns Returns the coordiantes of the given point
   in the @em reference space
 */

//////////////////////
// Inline Functions //
//////////////////////

inline double Jacobian::det(const fullVector<double>& UVW) const{
  return element->getJacobianDeterminant(UVW(0), 
					 UVW(1),
					 UVW(2));
}

inline double Jacobian::det(void) const{
  return 42;
}


#endif
