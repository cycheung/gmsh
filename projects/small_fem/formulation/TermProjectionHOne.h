#ifndef _TERMPROJECTIONHONE_H_
#define _TERMPROJECTIONHONE_H_

#include <vector>

#include "Jacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @interface TermProjectionHOne
   @brief Term for @f$ H^1 @f$ Projection

   Term for @f$ H^1 @f$ Projection
 */

class TermProjectionHOne: public Term{
 private:
  // Function to Project //
  double (*f)(fullVector<double>& xyz);

  // Integration Points //
  const fullVector<double>* gW;
  const fullMatrix<double>* gC;
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
  TermProjectionHOne(const Jacobian& jac,
                     const Basis& basis,
                     const fullVector<double>& integrationWeights,
                     const fullMatrix<double>& integrationPoints,
                     double (*f)(fullVector<double>& xyz));

  virtual ~TermProjectionHOne(void);

 private:
  void clean(void);
  void buildEMap(void);

  void computeC(void);
  void computeB(void);
  void computeA(void);
};

#endif
