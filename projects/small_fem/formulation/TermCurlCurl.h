#ifndef _TERMCURLCURL_H_
#define _TERMCURLCURL_H_

#include "GroupOfJacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermCurlCurl
   @brief Term of the @c Curl @c Curl type

   Term of the @c Curl @c Curl type
 */

class TermCurlCurl: public Term{
 private:
  typedef const fullMatrix<double>& (Basis::*bFunction)(unsigned int s)const;

 public:
  TermCurlCurl(const GroupOfJacobian& goj,
               const Basis& basis,
               const fullVector<double>& integrationWeights);

  virtual ~TermCurlCurl(void);

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
