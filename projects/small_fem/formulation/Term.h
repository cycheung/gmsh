#ifndef _TERM_H_
#define _TERM_H_

#include <map>

#include "fullMatrix.h"
#include "GroupOfDof.h"

/**
   @interface Term
   @brief Interface of helper methods for computing Finit Element Terms

   Interface of helper methods for computing Finit Element Terms
 */

class Term{
 protected:
  unsigned int nFunction;
  unsigned int nOrientation;
  std::vector<unsigned int>* orientationStat;

  fullMatrix<double>** aM;

 public:
  Term(void);
  virtual ~Term(void);

  double getTerm(unsigned int dofI,
                 unsigned int dofJ,
                 unsigned int elementId) const;
};

#endif
