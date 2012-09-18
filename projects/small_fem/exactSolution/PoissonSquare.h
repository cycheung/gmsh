#ifndef _POISSONSQUARE_H_
#define _POISSONSQUARE_H_

#include "ExactSolution.h"

/**
   @class PoissonSquare
   @brief Exact solution of Poisson equation over a square

   Exact soltuion of Poisson equation over a square.

   @f[
   \left\{
   \begin{array}{lllr}
   \Delta f(x, y) & = & 1 & \text{Over } \Omega\\
   f(x, y)        & = & 0 & \text{Over } \partial\Omega
   \end{array}
   \ritgh.
   @f]

   With @f$ \Omega = \Big{(x, y)\in [-1, 1]~\times~[-1, 1]\Big} @f$

   @note
   PoissonSquare is an ExactSolution
*/


class PoissonSquare: public ExactSolution{
 private:
  static double       pi;
  static unsigned int max;

 public:
  PoissonSquare(const GroupOfElement& goe);
  virtual ~PoissonSquare(void);

 private:
  virtual double fScalar(double x, double y, double z);
};


/**
   @fn PoissonSquare::PoissonSquare
   @param goe The Meshed square domain

   Instanciates a new PoissonSquare and 
   compute solution
   **

   @fn PoissonSquare::~PoissonSquare
   Deletes this PoissonSquare   
   **
*/

#endif
