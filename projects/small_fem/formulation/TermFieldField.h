#ifndef _TERMFIELDFIELD_H_
#define _TERMFIELDFIELD_H_

#include "GroupOfJacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermFieldField
   @brief A Term of the Field Field type

   A Term of the Field Field type
 */

class TermFieldField: public Term{
 public:
  TermFieldField(const GroupOfJacobian& goj,
                 const Basis& basis,
                 const fullVector<double>& integrationWeights);

  virtual ~TermFieldField(void);

 private:
  void computeC(const Basis& basis,
                const fullVector<double>& gW,
                fullMatrix<double>**& cM);

  void computeB(const GroupOfJacobian& goj,
                unsigned int nG,
                fullMatrix<double>**& bM);
};

#endif
