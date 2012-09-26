#ifndef _POISSONCIRCLE_H_
#define _POISSONCIRCLE_H_

#include "ExactSolution.h"

/**
   @class PoissonCircle
   @brief Exact solution of Poisson equation over a circle

   Exact soltuion of Poisson equation over a circle.

   @f[
   \left\{
   \begin{array}{lllr}
   \Delta f(x, y) & = & 1 & \text{Over } \Omega\\
   f(x, y)        & = & 0 & \text{Over } \partial\Omega
   \end{array}
   \ritgh.
   @f]

   With @f$ \Omega = \Big{(x, y) \Big| x^2 + y^2 \in [-1, 1]~\times~[-1, 1]\Big} @f$

   @note
   PoissonCircle is an ExactSolution
*/


class PoissonCircle: public ExactSolution{
 public:
  PoissonCircle(const GroupOfElement& goe);
  virtual ~PoissonCircle(void);

 private:
  virtual double fScalar(double x, double y, double z);
};


/**
   @fn PoissonCircle::PoissonCircle
   @param goe The Meshed circle domain

   Instanciates a new PoissonCircle and 
   compute solution
   **

   @fn PoissonCircle::~PoissonCircle
   Deletes this PoissonCircle   
   **
*/

#endif
