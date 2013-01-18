#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include <string>
#include <vector>

#include "System.h"
#include "EigenSystem.h"
#include "Writer.h"
#include "fullMatrix.h"

#include "Mesh.h"
#include "GroupOfElement.h"
#include "DofManager.h"
#include "FunctionSpace.h"

/**
   @class Interpolator
   @brief Interpolating methods for functions and FEM Solutions

   This class allows the interpolation of
   functions and FEM Solutions.@n

   The interpolation is done on the @em nodes of a given GroupOfElement.@n
   This GroupOfElement is defined when the Interpolator is Instanciated

   This class can also @em write the interpolation into a file.@n
   The file format is defined by a Writer.

   @todo
   Allow multiple Basis.
 */

class Interpolator{
 private:
  bool ownSol;
  const fullVector<double>* sol;

  const Mesh*          mesh;
  const DofManager*    dofM;
  const FunctionSpace* fs;

  double             (*fScalar)(fullVector<double>& xyz);
  fullVector<double> (*fVector)(fullVector<double>& xyz);

  bool                              scalar;
  std::vector<double>*              nodalScalarValue;
  std::vector<fullVector<double> >* nodalVectorValue;
  const GroupOfElement*             visuDomain;

 public:
   Interpolator(const System& system);
   Interpolator(const System& system,
		const GroupOfElement& visu);

   Interpolator(const EigenSystem& system,
		unsigned int eigenNumber);
   Interpolator(const EigenSystem& system,
		unsigned int eigenNumber,
		const GroupOfElement& visu);

   Interpolator(double (*f)(fullVector<double>& xyz),
	    const GroupOfElement& visu);
   Interpolator(fullVector<double> (*f)(fullVector<double>& xyz),
	    const GroupOfElement& visu);

  ~Interpolator(void);

  void write(const std::string name, Writer& writer) const;
  bool isScalar(void) const;

  std::vector<double>&              getNodalScalarValue(void) const;
  std::vector<fullVector<double> >& getNodalVectorValue(void) const;

 private:
  void initSystem(const System& system);
  void initSystem(const EigenSystem& system,
		  unsigned int eigenNumber);

  void interpolate(void);
  void interpolateOnVisu(void);
  void evaluateF(void);

  static const fullVector<double>*
    getSol(const std::vector<std::vector<std::complex<double> > >& eVector,
	   unsigned int eigenNumber);
};


/**
   @fn Interpolator::Interpolator(const System&)
   @param system A System

   Instanciate a new Interpolator,
   based on the given System%'s @em Solution

   @note
   The interpolation will be done on the  @em support of the
   System%'s Formulation%'s FunctionSpace.
   **

   @fn Interpolator::Interpolator(const System&, const GroupOfElement&)
   @param system A System
   @param visu The GroupOfElement to use for interpolation

   Instanciate a new Interpolator,
   based on the given System%'s @em Soltuion

   @note
   The interpolation will be done on the given GroupOfElement
   **

   @fn Interpolator::Interpolator(const EigenSystem&, unsigned int)
   @param system An EigenSystem
   @param eigenNumber A natural number

   Instanciate a new Interpolator,
   based on the given EigenSystem%'s Eigenvector
   number @c eigenNumber

   @note
   The interpolation will be done on the  @em support of the
   EigenSystem%'s EigenFormulation%'s FunctionSpace.
   **

   @fn Interpolator::Interpolator(const EigenSystem&, unsigned int, const GroupOfElement&)
   @param system A System
   @param eigenNumber a natural number
   @param visu The GroupOfElement to use for interpolation

   Instanciate a new Interpolator,
   based on the given EigenSystem%'s Eigenvector
   number @c eigenNumber

   @note
   The interpolation will be done on the given GroupOfElement
   **

   @fn Interpolator::Interpolator(double (*f)(fullVector<double>& xyz), const GroupOfElement&)
   @param f A @em scalar Function
   @param visu The GroupOfElement to use for interpolation

   Instanciate a new Interpolator,
   based on the given Funtion

   @note
   The interpolation will be done on the given GroupOfElement
   **

   @fn Interpolator::Interpolator(fullVector<double> (*f)(fullVector<double>& xyz), const GroupOfElement&)
   @param f A @em vectorial Function
   @param visu The GroupOfElement to use for interpolation

   Instanciate a new Interpolator,
   based on the given Funtion

   @note
   The interpolation will be done on the given GroupOfElement
   **

   @fn Interpolator::~Interpolator
   Deletes this Interpolator
   **

   @fn Interpolator::write
   @param name The file (@em without extension) where the
   interpolated values will be written
   @param writer The Writer to use

   Writes the interpolated values into the given file@n

   The file format is given by the Writer
   **

   @fn Interpolator::isScalar
   @return Returns:
   @li @c true, if the interpolated values are @em scalar
   @li @c false, otherwise
   **

   @fn Interpolator::getNodalScalarValue
   @return Returns the Interpolated @em scalar values
   @note If Interpolator::isScalar() is @c false
   an Exception will be thrown
   **

   @fn Interpolator::getNodalVectorValue
   @return Returns the Interpolated @em vectorial values
   @note If Interpolator::isScalar() is @c true
   an Exception will be thrown
 */

//////////////////////
// Inline Functions //
//////////////////////

inline bool Interpolator::isScalar(void) const{
  return scalar;
}

#endif
