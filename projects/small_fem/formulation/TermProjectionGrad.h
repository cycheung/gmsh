#ifndef _TERMPROJECTIONGRAD_H_
#define _TERMPROJECTIONGRAD_H_

#include "GroupOfJacobian.h"
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
  typedef const fullMatrix<double>& (Basis::*bFunction)(unsigned int s)const;

 public:
  TermProjectionGrad(const GroupOfJacobian& goj,
                     const Basis& basis,
                     const fullVector<double>& integrationWeights,
                     const fullMatrix<double>& integrationPoints,
                     fullVector<double> (*f)(fullVector<double>& xyz));

  virtual ~TermProjectionGrad(void);

 private:
  void computeC(const Basis& basis,
                const bFunction& getFunction,
                const fullVector<double>& gW,
                fullMatrix<double>**& cM);

  void computeB(const GroupOfJacobian& goj,
                const fullMatrix<double>& gC,
                fullVector<double> (*f)(fullVector<double>& xyz),
                fullMatrix<double>**& bM);
};

#endif
