#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "MElement.h"
#include "DofManager.h"
#include "FunctionSpace.h"
#include "fullMatrix.h"
#include "System.h"

/**
   @class Integrator
   @brief Integrate a FEM solution
   
   Compute the integral of a @em solved
   FEM System.
 */

class Integrator{
 private:
  // FEM Solution
  const fullVector<double>* sol;
  const DofManager*         dofM;
  const FunctionSpace*      fs;

  // Geometry
  const std::vector<const MElement*>*  element;
  unsigned int                        nElement;

  // Gaussian Quadrature Data
  fullMatrix<double>* gC;
  fullVector<double>* gW;

 public:
   Integrator(const System& system);
  ~Integrator(void);

  double integrate(void) const;

 private:
  double integrate(double (*law)(MElement& element, 
				 const FunctionSpace& fs, 
				 std::vector<double>& coef,
				 fullMatrix<double>& gC,
				 fullVector<double>& gW)) const;

  static double zero(MElement& element, 
		     const FunctionSpace& fs, 
		     std::vector<double>& coef,
		     fullMatrix<double>& gC,
		     fullVector<double>& gW);
  
  static double one(MElement& element,
		    const FunctionSpace& fs, 
		    std::vector<double>& coef,
		    fullMatrix<double>& gC,
		    fullVector<double>& gW);
  
  static double two(MElement& element, 
		    const FunctionSpace& fs, 
		    std::vector<double>& coef,
		    fullMatrix<double>& gC,
		    fullVector<double>& gW);
  
  static double three(MElement& element,
		      const FunctionSpace& fs, 
		      std::vector<double>& coef,
		      fullMatrix<double>& gC,
		      fullVector<double>& gW);
};


/**
   @fn Integrator::Integrator
   @param system A Finite Element System

   Instanciates a new Integrator that will be able
   to integrate the given System%'s solution
   **

   @fn Integrator::~Integrator
   
   Deletes this Integrator
   **

   @fn Integrator::integrate
   @return Returns the integrated solution
 */

#endif
