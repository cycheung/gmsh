#ifndef _FORMULATIONEIGENFREQUENCY_H_
#define _FORMULATIONEIGENFREQUENCY_H_

#include <vector>

#include "FunctionSpaceEdge.h"
#include "EigenFormulation.h"

/**
   @class FormulationEigenFrequency
   @brief Formulation for the Eigenfrequencies Problem

   Formulation for the Eigenfrequencies Problem
 */

class FormulationEigenFrequency: public EigenFormulation{
 private:
  // Physical Values //
  static const double mu;
  static const double eps;

  // Gaussian Quadrature Data (Term One) //
  int G1;
  fullMatrix<double>* gC1;
  fullVector<double>* gW1;

  // Gaussian Quadrature Data (Term Two) //
  int G2;
  fullMatrix<double>* gC2;
  fullVector<double>* gW2;

  // Function Space //
  FunctionSpaceEdge* fspace;

 public:
  FormulationEigenFrequency(const GroupOfElement& goe,
			    unsigned int order);

  virtual ~FormulationEigenFrequency(void);
  
  virtual double weakA(int dofI, int dofJ,
		       const GroupOfDof& god) const;

  virtual double weakB(int dofI, int dofJ,
		       const GroupOfDof& god) const;
  
  virtual double rhs(int equationI, 
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationEigenFrequency::FormulationEigenFrequency
   Instanciates a new EigenFormulation for the
   Eigenfrequencies Problem
   **

   @fn FormulationEigenFrequency::~FormulationEigenFrequency
   Deletes this FormualtionEigenFrequency
*/

//////////////////////
// Inline Functions //
//////////////////////

inline double FormulationEigenFrequency::rhs(int equationI, 
					     const GroupOfDof& god) const{
  return 0;
}

inline const FunctionSpace& FormulationEigenFrequency::fs(void) const{
  return *fspace;
}

#endif
