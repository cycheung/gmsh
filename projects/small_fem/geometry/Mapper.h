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

  static void hCurl(const fullMatrix<double>& hCurlUVW,
                    unsigned int              row,
                    unsigned int              col,
                    const fullMatrix<double>& invJac,
                    fullVector<double>&       hCurlXYZ);

  static void hDiv(const fullMatrix<double>& hDivUVW,
                   unsigned int              row,
                   unsigned int              col,
                   const fullMatrix<double>& jac,
                   double                    det,
                   fullVector<double>&       hDivXYZ);
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
 */

/*
   @fn Mapper::grad(const fullVector<double>& gradUVW, const fullMatrix<double>& invJac)
   @param gradUVW A gradient in the @em reference space
   @param invJac The Invert Jacobian Matrix evaluated at @c UVW
   @returns Returns the given gradient in the
   @em physical space
   **

   @fn Mapper::curl(const fullVector<double>& curlUVW, const fullMatrix<double>& jac, double invDet);
   @param curlUVW A curl in the @em reference space
   @param jac The Jacobian Matrix evaluated at @c UVW
   @param invDet The Invert of the Jacobian Matrix Determinant
   evaluated at @c UVW
   @returns Returns the given curl in the
   @em physical space
 */

#endif
