#ifndef _MAPPER_H_
#define _MAPPER_H_

#include "fullMatrix.h"

/**
   @class Mapper
   @brief Set of methods for mapping

   This class provides a set of class methods
   for handling mapping between physical and reference spaces.

   Because this class got only class method, it doesn't need to be instanciated.

   The pysical space is defined by the X, Y and Z coordinates.

   The reference space is defined by the U, V and W coordinates.
*/

class Mapper{
 public:
   Mapper(void);
  ~Mapper(void);

  static void hCurl(const fullMatrix<double>& hCurlUVW,
                    size_t                    row,
                    size_t                    col,
                    const fullMatrix<double>& invJac,
                    fullVector<double>&       hCurlXYZ);

  static void hDiv(const fullMatrix<double>& hDivUVW,
                   size_t                    row,
                   size_t                    col,
                   const fullMatrix<double>& jac,
                   double                    det,
                   fullVector<double>&       hDivXYZ);
};

/**
   @fn Mapper::Mapper
   Instanciates a new Mapper

   Mapper got only static methods, so it is not requiered to instanciate it
   **

   @fn Mapper::~Mapper
   Deletes this Mapper
   **

   @fn Mapper::hCurl
   @param hCurlUVW A set of @f$H(\mathrm{\mathbf{curl}})@f$ fields
   (each line is a field
   and each 3 columns (for the 3 coordinates) is a vector of this field)
   @param row A row index of hCurlUVW
   @param col A column triplet index of hCurlUVW
   @param invJac A invert jacobian matrix evaluated at UVW
   @param hCurlXYZ An allocated 3D vector that will contain the mapped vector

   Fills hCurlXYZ with the mapping (from the UVW to the XYZ space)
   of a vector given by
   [hCurlUVW(row, col * 3 + 0);
    hCurlUVW(row, col * 3 + 1);
    hCurlUVW(row, col * 3 + 2)]
   **

   @fn Mapper::hDiv
   @param hDivUVW A set of @f$H(\mathrm{div})@f$ fields
   (each line is a field
   and each 3 columns (for the 3 coordinates) is a vector of this field)
   @param row A row index of hDivUVW
   @param col A column triplet index of hDivUVW
   @param jac A jacobian matrix evaluated at UVW
   @param det The determinent of the given matrix
   @param hDivXYZ An allocated 3D vector that will contain the mapped vector

   Fills hDivXYZ with the mapping (from the UVW to the XYZ space)
   of a vector given by
   [hDivUVW(row, col * 3 + 0);
    hDivUVW(row, col * 3 + 1);
    hDivUVW(row, col * 3 + 2)]
 */

#endif
