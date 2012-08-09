#ifndef _FORMULATIONLAPLACE_H_
#define _FORMULATIONLAPLACE_H_

#include <vector>
#include "Formulation.h"
#include "Mesh.h"
#include "Polynomial.h"

/**
   @class FormulationLaplace
   @brief Formulation for the Laplace problem

   Formulation for the @em Laplace problem.

   @todo
   Remove const_cast by correcting MElement constness@n
   Allow Hybrid Mesh
 */

class FormulationLaplace: public Formulation{
 private:
  // Gaussian Quadrature Data //
  int G;
  fullMatrix<double>* gC;
  fullVector<double>* gW;

  // Function Space //
  FunctionSpace* fspace;

  // Grad Field //
  unsigned int             basisSize;
  std::vector<Polynomial>* gradBasis;

 public:
  FormulationLaplace(const GroupOfElement& goe, Mesh& mesh);

  virtual ~FormulationLaplace(void);

  virtual double weak(const int nodeI, const int nodeJ, 
		      const GroupOfDof& god) const;

  virtual double rhs(const int equationI,
		     const GroupOfDof& god) const;

  virtual FunctionSpace& fs(void) const;
};

/**
   @fn FormulationLaplace::FormulationLaplace
   @return Returns a new FormulationLaplace
 
   @fn FormulationLaplace::~FormulationLaplace
   @return Deletes this FormulationLaplace
*/

//////////////////////
// Inline Functions //
//////////////////////

inline double FormulationLaplace::rhs(const int equationI,
				      const GroupOfDof& god) const{
  return 0;
}

inline FunctionSpace& FormulationLaplace::fs(void) const{
  return *fspace;
}

#endif
