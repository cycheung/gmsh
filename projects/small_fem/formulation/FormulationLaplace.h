#ifndef _FORMULATIONLAPLACE_H_
#define _FORMULATIONLAPLACE_H_

#include <map>
#include <vector>

#include "FunctionSpaceScalar.h"
#include "Jacobian.h"
#include "Formulation.h"

/**
   @class FormulationLaplace
   @brief Formulation for the Laplace problem

   Formulation for the @em Laplace problem.

   @todo
   Remove ALL const_cast%S by correcting MElement constness@n
   Allow Hybrid Mesh
 */

class FormulationLaplace: public Formulation{
 private:
  // Gaussian Quadrature Data //
  int G;
  fullMatrix<double>* gC;
  fullVector<double>* gW;

  // Function Space & Basis//
  FunctionSpaceScalar* fspace;
  Basis*               basis;

  // 'Fast' Assembly //
  unsigned int    nOrientation;
  unsigned int    nFunction;
  unsigned int    nElement;
  GroupOfElement* goe;
  Jacobian*       jac;

  std::map<const MElement*, std::pair<unsigned int, unsigned int> >* eMap;
  std::vector<unsigned int>* orientationStat;

  fullMatrix<double>** cM;
  fullMatrix<double>** bM;
  fullMatrix<double>** aM;

 public:
  FormulationLaplace(GroupOfElement& goe, unsigned int order);

  virtual ~FormulationLaplace(void);

  virtual double weak(int nodeI, int nodeJ,
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
   @fn FormulationLaplace::FormulationLaplace
   @param goe A GroupOfElement
   @param order A natural number

   Instantiates a new FormulationLaplace of the given order@n

   The given GroupOfElement will be used as the
   geomtrical @em domain
   **

   @fn FormulationLaplace::~FormulationLaplace
   Deletes this FormulationLaplace
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

inline double FormulationLaplace::rhs(int equationI,
				      const GroupOfDof& god) const{
  return 0;
}

inline const FunctionSpace& FormulationLaplace::fs(void) const{
  return *fspace;
}

#endif
