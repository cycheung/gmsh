#ifndef _TERM_H_
#define _TERM_H_

#include "fullMatrix.h"

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

  mutable bool         once;
  mutable unsigned int lastId;
  mutable unsigned int lastI;
  mutable unsigned int lastCtr;

 public:
  Term(void);
  virtual ~Term(void);

  double getTerm(unsigned int dofI,
                 unsigned int dofJ,
                 unsigned int elementId) const;

 private:
  double getTermOutCache(unsigned int dofI,
                         unsigned int dofJ,
                         unsigned int elementId) const;

 protected:
  void computeA(fullMatrix<double>**& bM,
                fullMatrix<double>**& cM);

  void clean(fullMatrix<double>**& bM,
             fullMatrix<double>**& cM);
};

/////////////////////
// Inline Function //
/////////////////////

inline double Term::getTerm(unsigned int dofI,
                            unsigned int dofJ,
                            unsigned int elementId) const{

  if(!once || elementId != lastId)
    // If Out Of Cache --> Fetch
    return getTermOutCache(dofI, dofJ, elementId);

  else
    // Else, rock baby yeah !
    return (*aM[lastI])(lastCtr, dofI * nFunction + dofJ);
}

#endif
