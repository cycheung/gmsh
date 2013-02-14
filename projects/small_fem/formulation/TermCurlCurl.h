#ifndef _TERMCURLCURL_H_
#define _TERMCURLCURL_H_

#include <vector>

#include "Jacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermCurlCurl
   @brief Term of the @c Curl @c Curl type

   Term of the @c Curl @c Curl type
 */

class TermCurlCurl: public Term{
 private:
  // Integration Points //
  const fullVector<double>* gW;
  unsigned int              nG;

  // Basis & Jacobians //
  const fullMatrix<double>** phi;
  const Jacobian* jac;

  // FE Matrix
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

 public:
  TermCurlCurl(const Jacobian& jac,
               const Basis& basis,
               const fullVector<double>& integrationWeights);

  virtual ~TermCurlCurl(void);

 private:
  void clean(void);

  void computeC(void);
  void computeB(void);
  void computeA(void);
};

#endif
