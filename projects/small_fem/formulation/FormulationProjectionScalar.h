#ifndef _FORMULATIONPROJECTIONSCALAR_H_
#define _FORMULATIONPROJECTIONSCALAR_H_

#include <map>

#include "FunctionSpaceScalar.h"
#include "fullMatrix.h"
#include "Jacobian.h"
#include "Formulation.h"

/**
   @class FormulationProjectionScalar
   @brief Formulation for the Projection of a Scalar Function problem

   Scalar Formulation for the @em L2 @em Projection problem.
 */

class FormulationProjectionScalar: public Formulation{
 private:
  // Gaussian Quadrature Data //
  int G;
  fullMatrix<double>* gC;
  fullVector<double>* gW;

  // Function Space & Basis //
  FunctionSpaceScalar* fspace;
  const Basis*         basis;

  // Function to Project //
  double (*f)(fullVector<double>& xyz);

  // 'Fast' Assembly //
  unsigned int    nOrientation;
  unsigned int    nFunction;
  unsigned int    nElement;
  GroupOfElement* goe;
  Jacobian*       jac;

  std::map<const MElement*, std::pair<unsigned int, unsigned int> >* eMap;
  std::vector<unsigned int>* orientationStat;

  fullMatrix<double>** lcM;
  fullMatrix<double>** lbM;
  fullMatrix<double>** laM;

  fullMatrix<double>** rcM;
  fullMatrix<double>** rbM;
  fullMatrix<double>** raM;

 public:
  FormulationProjectionScalar(double (*f)(fullVector<double>& xyz),
			      FunctionSpaceScalar& fs);

  virtual ~FormulationProjectionScalar(void);

  virtual double weak(int dofI, int dofJ,
		      const GroupOfDof& god) const;

  virtual double rhs(int equationI,
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;

 private:
  void computeC(void);
  void computeB(void);
  void computeA(void);

  void deleteCB(void);
};

/**
   @fn FormulationProjectionScalar::FormulationProjectionScalar
   @param f The function to project
   @param fs A FunctionSpaceNode

   Instantiates a new FormulationProjectionScalar to project
   the given function@n

   FormulationProjectionScalar will use the given FunctionSpace
   for the projection

   @warning The given FunctionSpace will be ask to
   @em PreEvaluate
   (@em see FunctionSpaceScalar::preEvaluateLocalFunctions()
   and similar)
   some Functions
   **

   @fn FormulationProjectionScalar::~FormulationProjectionScalar
   Deletes the this FormulationProjectionScalar
   **
*/


//////////////////////
// Inline Functions //
//////////////////////

inline const FunctionSpace& FormulationProjectionScalar::fs(void) const{
  return *fspace;
}

#endif
