#ifndef _TERMFIELDFIELD_H_
#define _TERMFIELDFIELD_H_

#include <vector>

#include "Jacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermFieldField
   @brief A Term of the Field Field type

   A Term of the Field Field type
 */

class TermFieldField: public Term{
 private:
  // Integration Points //
  const fullVector<double>* gW;
  unsigned int              nG;

  // Basis & Jacobians //
  const Basis*    basis;
  const Jacobian* jac;

  // FE Matrix
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

 public:
  TermFieldField(const Jacobian& jac,
                 const Basis& basis,
                 const fullVector<double>& integrationWeights);

  virtual ~TermFieldField(void);

 private:
  void clean(void);

  void computeC(void);
  void computeB(void);
  void computeA(void);
};

#endif
