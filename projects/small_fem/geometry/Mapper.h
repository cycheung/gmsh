#ifndef _MAPPER_H_
#define _MAPPER_H_

#include "fullMatrix.h"

/**
   @class Mapper
   @brief Set of methods for mapping

   This class provides a set of @em class @em methods
   for handling mapping between physical 
   and reference spaces.@n

   @note
   Because this class got @em only @em class @em method,
   it @em doesn't need to be instanciated.

   The @em pysical space is defined by the
   @c X, @c Y and @c Z coordinates.@n

   The @em reference space is defined by the
   @c U, @c V and @c W coordinates.@n
*/

class Mapper{
 public:
   Mapper(void);
  ~Mapper(void);

  static fullVector<double> map(const fullVector<double>& UVW, 
				const fullVector<double>& originXYZ,
				const fullMatrix<double>& jac);

  static fullVector<double> invMap(const fullVector<double>& XYZ, 
				   const fullVector<double>& originXYZ,
				   const fullMatrix<double>& invJac);

  static fullVector<double> grad(const fullVector<double>& gradUVW, 
				 const fullMatrix<double>& invJac);
};

/**
   @fn Mapper::Mapper
   Instanciates a new Mapper
   @note Mapper got @em only @em class
   methods (functions), so it is not requiered
   to instanciate a Mapper
   **

   @fn Mapper::~Mapper
   Deletes this Mapper
   **

   @fn Mapper::map
   @param UVW A @c 3D Vector with the coordinate 
   of a point in the @em reference space
   @param originXYZ A @c 3D Vector with the coordinate 
   (in the Physical Space) of a the first vertex of 
   the mapped element
   @param jac The Jacobian Matrix evaluated at @c UVW 
   @returns Returns the coordiantes of the given point
   in the @em physical space

   @fn Mapper::invMap
   @param XYZ A @c 3D Vector with the coordinate 
   of a point in the @em physical space
   @param originXYZ A @c 3D Vector with the coordinate 
   (in the Physical Space) of a the first vertex of 
   the mapped element
   @param invJac The Invert Jacobian Matrix evaluated at @c XYZ 
   @returns Returns the coordiantes of the given point
   in the @em reference space

   @fn Mapper::grad
   @param gradUVW A gradient in the @em reference space
   @param invJac The Invert Jacobian Matrix evaluated at @c UVW 
   @returns Returns the given gradient in the 
   @em physical space
 */

#endif
