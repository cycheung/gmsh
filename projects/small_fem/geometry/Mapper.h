#ifndef _MAPPER_H_
#define _MAPPER_H_

#include "fullMatrix.h"
#include "MElement.h"

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

   @todo
   Implement a cache ?@n
   Change MElement with const qualifier, so we can use
   'const MElement' in prototypes
*/

class Mapper{
 public:
   Mapper(void);
  ~Mapper(void);

  static double det(double u, double v, double w,
		    MElement& element);

  static double det(const fullVector<double>& UVW, 
		    MElement& element);

  static fullVector<double> map(const fullVector<double>& UVW, 
				MElement& element);

  static fullVector<double> grad(const fullVector<double>& gradUVW, 
				 MElement& element);

  static fullVector<double> invMap(const fullVector<double>& XYZ, 
				   MElement& element);
};

/**
   @fn Mapper::Mapper
   @return Returns a new Mapper
   @note Mapper has @em only @em static
   methods (functions), so it is not requiered
   to instanciate a Mapper.

   @fn Mapper::~Mapper
   @return Deletes this Mapper

   @fn Mapper::det(const fullVector<double>&) 
   @param UVW A @c 3D Vector with the coordinate
   of a point in the @em reference space
   @return Returns the determinant of the 
   jacobian matrix at the given point

   @fn Mapper::map(const fullVector<double>&) 
   @param UVW A @c 3D Vector with the coordinate 
   of a point in the @em reference space
   @returns Returns the coordiantes of the given point
   in the @em physical space

   @fn Mapper::grad(const fullVector<double>&) 
   @param gradUVW A gradient in the @em reference space
   @returns Returns the given gradient in the 
   @em physical space

   @fn Mapper::invMap(const fullVector<double>&) 
   @param XYZ A @c 3D Vector with the coordinate 
   of a point in the @em physical space
   @returns Returns the coordiantes of the given point
   in the @em reference space
 */

//////////////////////
// Inline Functions //
//////////////////////

inline double Mapper::det(const fullVector<double>& UVW, MElement& element){
  return element.getJacobianDeterminant(UVW(0), 
					UVW(1),
					UVW(2));
}

inline double Mapper::det(double u, double v, double w, MElement& element){
  return element.getJacobianDeterminant(u, v, w); 
}


#endif
