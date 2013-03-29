#include <cmath>
#include "Mapper.h"

Mapper::Mapper(void){
}

Mapper::~Mapper(void){
}

//   WARNING                  //
// Jacobians are transposed ! //

void Mapper::hCurl(const fullMatrix<double>& hCurlUVW,
                   unsigned int              row,
                   unsigned int              col,
                   const fullMatrix<double>& invJac,
                   fullVector<double>&       hCurlXYZ){

  hCurlXYZ(0) =
    invJac(0, 0) * hCurlUVW(row, col * 3)     +
    invJac(0, 1) * hCurlUVW(row, col * 3 + 1) +
    invJac(0, 2) * hCurlUVW(row, col * 3 + 2);

  hCurlXYZ(1) =
    invJac(1, 0) * hCurlUVW(row, col * 3)     +
    invJac(1, 1) * hCurlUVW(row, col * 3 + 1) +
    invJac(1, 2) * hCurlUVW(row, col * 3 + 2);

  hCurlXYZ(2) =
    invJac(2, 0) * hCurlUVW(row, col * 3)     +
    invJac(2, 1) * hCurlUVW(row, col * 3 + 1) +
    invJac(2, 2) * hCurlUVW(row, col * 3 + 2);
}

void Mapper::hDiv(const fullMatrix<double>& hDivUVW,
                  unsigned int              row,
                  unsigned int              col,
                  const fullMatrix<double>& jac,
                  double                    det,
                  fullVector<double>&       hDivXYZ){

  hDivXYZ(0) =
    jac(0, 0) * hDivUVW(row, col * 3)     +
    jac(1, 0) * hDivUVW(row, col * 3 + 1) +
    jac(2, 0) * hDivUVW(row, col * 3 + 2);

  hDivXYZ(1) =
    jac(0, 1) * hDivUVW(row, col * 3)     +
    jac(1, 1) * hDivUVW(row, col * 3 + 1) +
    jac(2, 1) * hDivUVW(row, col * 3 + 2);

  hDivXYZ(2) =
    jac(0, 2) * hDivUVW(row, col * 3)     +
    jac(1, 2) * hDivUVW(row, col * 3 + 1) +
    jac(2, 2) * hDivUVW(row, col * 3 + 2);

  hDivXYZ(0) /= fabs(det);
  hDivXYZ(1) /= fabs(det);
  hDivXYZ(2) /= fabs(det);
}
