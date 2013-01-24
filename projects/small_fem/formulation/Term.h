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
  std::map<const MElement*, std::pair<unsigned int, unsigned int> >* eMap;
  unsigned int nFunction;

  fullMatrix<double>** aM;

 public:
  Term(void);
  virtual ~Term(void);

  double getTerm(unsigned int i,
                 unsigned int j,
                 const GroupOfDof& god) const;
};

//////////////////////
// Inline Functions //
//////////////////////

inline double Term::getTerm(unsigned int i,
                            unsigned int j,
                            const GroupOfDof& god) const{

  std::map<const MElement*, std::pair<unsigned int, unsigned int> >::iterator
    index = eMap->find(&god.getGeoElement());

  return (*aM[index->second.first])
    (index->second.second, i * nFunction + j);
}

#endif
