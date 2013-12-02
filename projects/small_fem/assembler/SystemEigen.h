#ifndef _SYSTEMEIGEN_H_
#define _SYSTEMEIGEN_H_

#include <complex>
#include "SystemTyped.h"
#include "petscmat.h"

/**
   @class SystemEigen
   @brief This class assembles an eigenvalue system

   This class assembles an eigenvalue system described by a Formulation.

   The eigenvalue problem can be generalized or not:
   @li An eigenvalue problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{I})\mathbf{x} = \mathbf{b}@f$
   @li A Generalized Eigenvalue problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{B})\mathbf{x} = \mathbf{b}@f$

   The Solver used is <a href="http://www.grycap.upv.es/slepc/">SLEPc</a>.
 */

class SystemEigen: public SystemTyped<std::complex<double> >{
 protected:
  bool general;

  Mat* A;
  Mat* B;

  PetscInt nEigenValues;
  fullVector<std::complex<double> >* eigenValue;
  std::vector<fullVector<std::complex<double> > >* eigenVector;

 public:
  SystemEigen(const FormulationTyped<std::complex<double> >& formulation);
  virtual ~SystemEigen(void);

  bool isGeneral(void) const;

  virtual size_t getNComputedSolution(void)                          const;
  virtual void   getSolution(fullVector<std::complex<double> >& sol,
                             size_t nSol)                            const;
  virtual void   getSolution(fullVector<std::complex<double> >& sol) const;

  void getEigenValues(fullVector<std::complex<double> >& eig)  const;

  void setNumberOfEigenValues(size_t nEigenValues);

  virtual void assemble(void);
  virtual void solve(void);

  virtual void addSolution(FEMSolution& feSol) const;
};


/**
   @fn SystemEigen::SystemEigen
   @param formulation A Formulation that
   gives the way to assemble the eigenvalue system

   Instantiates a new SystemEigen
   ***

   @fn SystemEigen::~SystemEigen
   Deletes this SystemEigen
   **

   @fn SystemEigen::isGeneral
   @return Returns:
   @li true, if the SystemEigen is a generalized one
   @li false otherwise
   **

   @fn SystemEigen::getEigenValues
   @param eig
   Allocate and populates eig with the eigenvalues of this SystemEigen
   **

   @fn SystemEigen::setNumberOfEigenValues
   @param nEigenValues A natural number

   Sets the number of eigenvalues computed by
   SystemEigen::solve() to the given number
*/

#endif
