#ifndef _TERMPROJECTIONGRAD_H_
#define _TERMPROJECTIONGRAD_H_

#include <vector>

#include "Jacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermProjectionGrad
   @brief Term of the Field (in physical space) Grad (in reference space) type

   Term of the Field (in physical space)
   Grad (in reference space) type
 */

class TermProjectionGrad: public Term{
 private:
  // Function to Project //
  fullVector<double> (*f)(fullVector<double>& xyz);

  // Integration Points //
  const fullVector<double>* gW;
  const fullMatrix<double>* gC;
  unsigned int              nG;

  // Basis & Jacobians //
  const fullMatrix<double>** phi;
  const Jacobian* jac;

  // FE Matrix
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

 public:
  TermProjectionGrad(const Jacobian& jac,
                     const Basis& basis,
                     const fullVector<double>& integrationWeights,
                     const fullMatrix<double>& integrationPoints,
                     fullVector<double> (*f)(fullVector<double>& xyz));

  virtual ~TermProjectionGrad(void);

 private:
  void clean(void);

  void computeC(void);
  void computeB(void);
  void computeA(void);
};

#endif
