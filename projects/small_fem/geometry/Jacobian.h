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
   Implement a cache ?
*/

class Jacobian{
 private:
  MElement* element;

 public:
   Jacobian(const MElement& element);
  ~Jacobian(void);

  double det(void) const;
  double det(const fullVector<double>& UV) const;

  fullVector<double> map(const fullVector<double>& UV) const;

  fullVector<double> grad(const fullVector<double>& gradUV) const;

  fullVector<double> invMap(const fullVector<double>& XY) const;
};

/**
   @fn Jacobian::Jacobian
   @param element The Element from which the Jacobian 
   will be computed
   @return Returns a new Jacobian

   @fn Jacobian::~Jacobian
   @return Deletes this Jacobian

   @fn Jacobian::det
   @return Returns the determinant of the 
   transformation jacobian matrix

   @fn Jacobian::grad(const fullVector<double>&) const
   @param gradUV A gradient in the @em reference space
   @returns Returns the given gradient in the 
   @em physical space

   @fn Jacobian::grad(const double, const double) const
   @param u The @c U coordinate 
   of a gradient in the @em reference (@c 2D) space
   @param v The @c V coordinate of the same gradient
   @returns Returns the given gradient in the 
   @em physical space

   @fn Jacobian::invMap(const fullVector<double>&) const
   @param XY A @c 2D Vector with the coordinates 
   of a point in the @em physical (@c 2D) space
   @returns Returns the coordiantes of the given point
   in the @em reference space

   @fn Jacobian::invMap(const double, const double) const
   @param x The @c X coordinate 
   of a point in the @em physical (@c 2D) space
   @param y The @c Y coordinate of the same point
   @returns Returns the coordiantes of the given point
   in the @em reference space

   @fn Jacobian::map(const fullVector<double>&) const
   @param UV A @c 2D Vector with the coordinates 
   of a point in the @em reference (@c 2D) space
   @returns Returns the coordiantes of the given point
   in the @em physical space

   @fn Jacobian::map(const double, const double) const
   @param u The @c U coordinate 
   of a point in the @em reference (@c 2D) space
   @param v The @c V coordinate of the same point
   @returns Returns the coordiantes of the given point
   in the @em physical space
 */

//////////////////////
// Inline Functions //
//////////////////////


#endif
