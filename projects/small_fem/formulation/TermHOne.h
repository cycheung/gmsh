#ifndef _TERMHONE_H_
#define _TERMHONE_H_

#include <vector>

#include "Jacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @interface TermHOne
   @brief Term for @f$ H^1 @f$ Terms

   Term for @f$ H^1 @f$ Terms
 */

class TermHOne: public Term{
 private:
  // Integration Points //
  const fullVector<double>* gW;
  unsigned int              nG;

  // Basis & Jacobians //
  const Basis*    basis;
  const Jacobian* jac;

  unsigned int               nOrientation;
  std::vector<unsigned int>* orientationStat;

  // FE Matrix
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

 public:
  TermHOne(const Jacobian& jac,
           const Basis& basis,
           const fullVector<double>& integrationWeights);

  virtual ~TermHOne(void);

 private:
  void clean(void);
  void buildEMap(void);

  void computeC(void);
  void computeB(void);
  void computeA(void);
};

#endif
