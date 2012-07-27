#ifndef _MAPPER_H_
#define _MAPPER_H_

#include "fullMatrix.h"

/**
   @class Mapper
   @brief Set of methods for mapping

   This class provides a set of @em static methods
   for handling mapping between physical 
   and reference spaces.@n

   The @em pysical space is defined by the
   @li @c X, @c Y and @c Z coordinates.@n

   The @em reference space is defined by the
   @li @c U, @c V and @c W coordinates.@n
*/

class Mapper{
 public:
   Mapper(void);
  ~Mapper(void);

  static fullVector<double> map(const fullVector<double>& UVW, 
				const fullMatrix<double>& jac);

  static fullVector<double> grad(const fullVector<double>& gradUVW, 
				 const fullMatrix<double>& invJac);

  static fullVector<double> invMap(const fullVector<double>& XYZ, 
				   const fullMatrix<double>& invJac);
};

/**
   @fn Mapper::Mapper
   @return Returns a new Mapper
   @note Mapper has @em only @em static
   methods (functions), so it is not requiered
   to instanciate a Mapper.

   @fn Mapper::~Mapper
   @return Deletes this Mapper

   @fn Mapper::map(const fullVector<double>&) 
   @param UVW A @c 3D Vector with the coordinate 
   of a point in the @em reference space
   @param jac The Jacobian Matrix evaluated at @c UVW 
   @returns Returns the coordiantes of the given point
   in the @em physical space

   @fn Mapper::grad(const fullVector<double>&) 
   @param gradUVW A gradient in the @em reference space
   @param invJac The Invert Jacobian Matrix evaluated at @c UVW 
   @returns Returns the given gradient in the 
   @em physical space

   @fn Mapper::invMap(const fullVector<double>&) 
   @param XYZ A @c 3D Vector with the coordinate 
   of a point in the @em physical space
   @param invJac The Invert Jacobian Matrix evaluated at @c UVW 
   @returns Returns the coordiantes of the given point
   in the @em reference space
 */

//////////////////////
// Inline Functions //
//////////////////////

inline fullVector<double> Mapper::map(const fullVector<double>& UVW,
				      const fullMatrix<double>& jac){
  fullVector<double> XYZ(3);
  jac.mult(UVW, XYZ);
  return XYZ;
}

inline fullVector<double> Mapper::grad(const fullVector<double>& gradUVW, 
				       const fullMatrix<double>& invJac){

  fullVector<double> gradXYZ(3);
  invJac.multWithATranspose(gradUVW, 1, 0, gradXYZ);
  return gradXYZ;
}

#endif
