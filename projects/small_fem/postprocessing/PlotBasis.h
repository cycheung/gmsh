#ifndef _PLOTBASIS_H_
#define _PLOTBASIS_H_

#include "Writer.h"

#include "GroupOfElement.h"
#include "fullMatrix.h"
#include "BasisLocal.h"

#include <string>
#include <vector>

/**
   @class PlotBasis
   @brief A Ploter for a Local Basis

   A PlotBasis allows the @em evaluation of @em all
   the functions of a Local Basis, on a given @em domain.@n

   A PlotBasis can write a file, that can be used to @em plot
   the Local Basis functions.@n

   The file format is given by a Writer.
 */

class PlotBasis{
 private:
  Writer*      writer;
  unsigned int nFunction;
  bool         isScalar;

  int N;
  int E;

  const std::vector<const MElement*>* element;
  const std::vector<MVertex*>*        node;

  std::vector<double>**              nodalScalarValue;
  std::vector<fullVector<double> >** nodalVectorValue;

 public:
  PlotBasis(const BasisLocal& basis,
	    const GroupOfElement& group,
	    Writer& writer);

  virtual ~PlotBasis(void);

  virtual void plot(const std::string name) const;

 private:
  void getGeometry(const GroupOfElement& group);
  void interpolateScalar(const BasisLocal& basis);
  void interpolateVector(const BasisLocal& basis);
};

/**
   @fn PlotBasis::PlotBasis
   @param basis The Local Basis functions to plot
   @param group A GroupOfElement,
   defining the @em geometrical @em domain of
   the Basis functions
   @param writer The Writer to use to write the output file

   Instantiates a new PlotBasis
   **

   @fn PlotBasis::~PlotBasis
   Deletes this PlotBasis
   **

   @fn PlotBasis::plot
   @param name The file to write the ploted Basis
   (@em without extension)

   Writes the ploted Basis functions in the
   given file
   **
 */

#endif
