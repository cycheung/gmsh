#ifndef _QUADRATURE_H_
#define _QUADRATURE_H_

#include "fullMatrix.h"

/**
   @class Quadrature
   @brief Quadrature points and weights

   An Integration Quadrature points and weights

   @note
   Since an order '1' Quadrangular Basis
   is bilinear (of order '2'), an order @c N
   Quadrature over a Quad will compute a quadrature
   such that a function of order @c 2 @c * @c N
   can be integrated
*/

class Quadrature{
 private:
  fullMatrix<double>* gC;
  fullVector<double>* gW;

 public:
   Quadrature(int elementType, int order,
              unsigned int multiplicity);

  ~Quadrature(void);

  const fullMatrix<double>& getPoints(void)  const;
  const fullVector<double>& getWeights(void) const;

 private:
  static void point(fullMatrix<double>& gC,
                    fullVector<double>& gW);
};

/**
   @fn Quadrature::Quadrature
   @param elementType An element type tag
   @param order An integer
   @param multiplicity A natural number

   Instantiates a new Quadrature
   over the requested element type and for
   function of order: @c order @c * @c multiplicty

   @note
   If @c order is zero or less, an order of
   @em one is assumed
   **

   @fn Quadrature::~Quadrature

   Deletes this Quadrature
   **

   @fn Quadrature::getPoints
   @return Returns a matrix with the integration points
   for this Quadrature

   @note The returned matrix got the following pattern:
   @li Each row is an integration points
   @li Each column is a dimension
   **

   @fn Quadrature::getWeights
   @return Returns a matrix with the integration weights
   for this Quadrature
*/

//////////////////////
// Inline Functions //
//////////////////////

inline const fullMatrix<double>&
Quadrature::getPoints(void) const{
  return *gC;
}

inline const fullVector<double>&
Quadrature::getWeights(void) const{
  return *gW;
}

#endif
