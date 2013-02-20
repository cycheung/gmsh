#ifndef _TERMGRADGRAD_H_
#define _TERMGRADGRAD_H_

#include "GroupOfJacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermGradGrad
   @brief Term of the @c Grad @c Grad type

   Term of the @c Grad @c Grad type
 */

class TermGradGrad: public Term{
 private:
  typedef const fullMatrix<double>& (Basis::*bFunction)(unsigned int s)const;

 public:
  TermGradGrad(const GroupOfJacobian& goj,
               const Basis& basis,
               const fullVector<double>& integrationWeights);

  virtual ~TermGradGrad(void);

 private:
  void computeC(const Basis& basis,
                const bFunction& getFunction,
                const fullVector<double>& gW,
                fullMatrix<double>**& cM);

  void computeB(const GroupOfJacobian& goj,
                unsigned int nG,
                fullMatrix<double>**& bM);
};

#endif
