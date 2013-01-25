#ifndef _TERMHDIV_H_
#define _TERMHDIV_H_

#include <vector>

#include "Jacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @interface TermHDiv
   @brief Term for @f$ H(\mathbf{\mathrm{div}}) @f$ Terms

   Term for @f$ H(\mathbf{\mathrm{div}}) @f$ Terms
 */

class TermHDiv: public Term{
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
  TermHDiv(const Jacobian& jac,
            const Basis& basis,
            const fullVector<double>& integrationWeights);

  virtual ~TermHDiv(void);

 private:
  void clean(void);
  void buildEMap(void);

  void computeC(void);
  void computeB(void);
  void computeA(void);
};

#endif
