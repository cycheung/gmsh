#ifndef _TERMPROJECTIONFIELD_H_
#define _TERMPROJECTIONFIELD_H_

#include <vector>

#include "Jacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermProjectionField
   @brief Term of the Field (in physical space) Field (in reference space) type

   Term of the Field (in physical space)
   Field (in reference space) type
 */

class TermProjectionField: public Term{
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

  // FE Matrix
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

 public:
  TermProjectionField(const Jacobian& jac,
                      const Basis& basis,
                      const fullVector<double>& integrationWeights,
                      const fullMatrix<double>& integrationPoints,
                      double (*f)(fullVector<double>& xyz));

  virtual ~TermProjectionField(void);

 private:
  void clean(void);

  void computeC(void);
  void computeB(void);
  void computeA(void);
};

#endif
