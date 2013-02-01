#ifndef _TERMGRADGRAD_H_
#define _TERMGRADGRAD_H_

#include <vector>

#include "Jacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermGradGrad
   @brief Term of the @c Grad @c Grad type

   Term of the @c Grad @c Grad type
 */

class TermGradGrad: public Term{
 private:
  // Integration Points //
  const fullVector<double>* gW;
  unsigned int              nG;

  // Basis & Jacobians //
  const fullMatrix<double>** phi;
  const Jacobian* jac;

  unsigned int               nOrientation;
  std::vector<unsigned int>* orientationStat;

  // FE Matrix
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

 public:
  TermGradGrad(const Jacobian& jac,
               const Basis& basis,
               const fullVector<double>& integrationWeights);

  virtual ~TermGradGrad(void);

 private:
  void clean(void);
  void buildEMap(void);

  void computeC(void);
  void computeB(void);
  void computeA(void);
};

#endif
