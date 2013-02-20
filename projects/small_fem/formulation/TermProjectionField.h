#ifndef _TERMPROJECTIONFIELD_H_
#define _TERMPROJECTIONFIELD_H_

#include "GroupOfJacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermProjectionField
   @brief Term of the Field (in physical space) Field (in reference space) type

   Term of the Field (in physical space)
   Field (in reference space) type
 */

class TermProjectionField: public Term{
 public:
  TermProjectionField(const GroupOfJacobian& goj,
                      const Basis& basis,
                      const fullVector<double>& integrationWeights,
                      const fullMatrix<double>& integrationPoints,
                      double (*f)(fullVector<double>& xyz));

  virtual ~TermProjectionField(void);

 private:
  void computeC(const Basis& basis,
                const fullVector<double>& gW,
                fullMatrix<double>**& cM);

  void computeB(const GroupOfJacobian& goj,
                const fullMatrix<double>& gC,
                double (*f)(fullVector<double>& xyz),
                fullMatrix<double>**& bM);
};

#endif
