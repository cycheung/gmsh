#ifndef _TERMPROJECTIONHCURL_H_
#define _TERMPROJECTIONHCURL_H_

#include <vector>

#include "Jacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @interface TermProjectionHCurl
   @brief Term for @f$ H(\mathbf{\mathrm{curl}}) @f$ Projections

   Term for @f$ H(\mathbf{\mathrm{curl}}) @f$ Projections
 */

class TermProjectionHCurl: public Term{
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

  unsigned int               nOrientation;
  std::vector<unsigned int>* orientationStat;

  // FE Matrix
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

 public:
  TermProjectionHCurl(const Jacobian& jac,
                      const Basis& basis,
                      const fullVector<double>& integrationWeights,
                      const fullMatrix<double>& integrationPoints,
                      fullVector<double> (*f)(fullVector<double>& xyz));

  virtual ~TermProjectionHCurl(void);

 private:
  void clean(void);
  void buildEMap(void);

  void computeC(void);
  void computeB(void);
  void computeA(void);
};

#endif
