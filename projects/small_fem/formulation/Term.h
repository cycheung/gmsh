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
  size_t nFunction;
  size_t nOrientation;
  const std::vector<size_t>* orientationStat;

  fullMatrix<double>** aM;

  mutable bool   once;
  mutable size_t lastId;
  mutable size_t lastI;
  mutable size_t lastCtr;

 public:
  virtual ~Term(void);

  double getTerm(size_t dofI,
                 size_t dofJ,
                 size_t elementId) const;

 private:
  double getTermOutCache(size_t dofI,
                         size_t dofJ,
                         size_t elementId) const;

 protected:
  Term(void);

  void allocA(size_t nFunction);

  void computeA(fullMatrix<double>**& bM,
                fullMatrix<double>**& cM);

  void clean(fullMatrix<double>**& bM,
             fullMatrix<double>**& cM);
};

/**
   @fn Term::~Term
   Deletes this Term
   **

   @fn Term::getTerm
   @param dofI A FunctionSpace function index
   @param dofJ A FunctionSpace function index
   @param elementId The ID of an Element
   @return Returns the finite element term associated to the given values
 */

/////////////////////
// Inline Function //
/////////////////////

inline double Term::getTerm(size_t dofI,
                            size_t dofJ,
                            size_t elementId) const{

  if(!once || elementId != lastId)
    // If Out Of Cache --> Fetch
    return getTermOutCache(dofI, dofJ, elementId);

  else
    // Else, rock baby yeah !
    return (*aM[lastI])(lastCtr, dofI * nFunction + dofJ);
}

#endif
