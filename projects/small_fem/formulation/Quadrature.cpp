#include "Quadrature.h"
#include "GaussIntegration.h"
#include "GmshDefines.h"
#include "Exception.h"

Quadrature::Quadrature(int elementType,
                       int order,
                       size_t multiplicity){
  // Alloc //
  gC = new fullMatrix<double>(1, 1);
  gW = new fullVector<double>(1);

  // In case of Nedelec Basis: use order 1 //
  if(order <= 0)
    order = 1;

  // True order:
  size_t trueOrder = order * multiplicity;

  // Get points and weigths //
  switch(elementType){
  case TYPE_PNT :
    point(*gC, *gW);
    break;

  case TYPE_LIN:
    gaussIntegration::get(TYPE_LIN, trueOrder, *gC, *gW);
    break;

  case TYPE_TRI:
    gaussIntegration::get(TYPE_TRI, trueOrder, *gC, *gW);
    break;

  case TYPE_QUA:
    gaussIntegration::get(TYPE_QUA, 2 * trueOrder, *gC, *gW);
    break;

  case TYPE_TET:
    gaussIntegration::get(TYPE_TET, trueOrder, *gC, *gW);
    break;

  case TYPE_HEX:
    gaussIntegration::get(TYPE_HEX, 3 * trueOrder, *gC, *gW);
    break;

  case TYPE_PRI:
    throw Exception("Quadrature on Prism not implemented");

  case TYPE_PYR:
    throw Exception("Quadrature on Pyramid not implemented");

  default:
    throw Exception("Quadrature: unknown element type (%d)", elementType);
  }
}

Quadrature::~Quadrature(void){
  delete gC;
  delete gW;
}

void Quadrature::point(fullMatrix<double>& gC,
                       fullVector<double>& gW){
  gW.resize(1);
  gC.resize(1, 1);

  gW(0)    = 1;
  gC(0, 0) = 0;
  gC(0, 1) = 0;
  gC(0, 2) = 0;
}
