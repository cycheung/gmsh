#ifndef _JACOBIAN_H_
#define _JACOBIAN_H_

#include <vector>
#include "Vector.h"
#include "fullMatrix.h"
#include "Node.h"

/**
   @class Jacobian
   @brief Transformations between
   physical and reference spaces

   This class handles the transformations
   between physical and reference spaces.@n

   The @em pysical space is defined by:
   @li @c X and @c Y coordinates in @c 2D
   @li @c X, @c Y and @c Z coordinates in @c 3D

   The @em reference space is defined by:
   @li @c U and @c V coordinates in @c 2D
   @li @c U, @c V and @c W coordinates in @c 3D

   @todo
   Use a real Matrix to handle Jacobian%s
*/

class Jacobian{
 private:
  int nNode;

  double* nodeX;
  double* nodeY;
  double* nodeZ;

  fullMatrix<double>* jac; // From Ref. to Phys. Space

  double dxdu;
  double dxdv;
  double dydu;
  double dydv;

  double detDxDu;

  double dudx;
  double dudy;
  double dvdx;
  double dvdy;

 public:
   Jacobian(const std::vector<Node*>& nodes);
  ~Jacobian(void);

  double det(void) const;

  Vector<double> grad(const Vector<double>& gradUV) const;

  Vector<double> invMap(const Vector<double>& XY) const;
  Vector<double> invMap(const double x, const double y) const;

  Vector<double> map(const Vector<double>& UV) const;
  Vector<double> map(const double u, const double v) const;

 private:
  void triJac(void);
};

/**
   @fn Jacobian::Jacobian
   @param nodes Node%s defining the geometry of the 
   @em physical element to transform (onto the 
   @em reference element)
   @return Returns a new Jacobian

   @fn Jacobian::~Jacobian
   @return Deletes this Jacobian

   @fn Jacobian::det
   @return Returns the determinant of the 
   transformation jacobian matrix

   @fn Jacobian::grad
   @param gradUV A gradient in the @em reference space
   @returns Returns the given gradient in the 
   @em physical space

   @fn Jacobian::invMap(const Vector<double>&) const
   @param XY A @c 2D Vector with the coordinates 
   of a point in the @em physical (@c 2D) space
   @returns Returns the coordiantes of the given point
   in the @em reference space

   @fn Jacobian::invMap(const double, const double) const
   @param x The @c X coordinate 
   of a point in the @em physical (@c 2D) space
   @param y The @c Y coordinate of the same point
   @returns Returns the coordiantes of the given point
   in the @em reference space

   @fn Jacobian::map(const Vector<double>&) const
   @param UV A @c 2D Vector with the coordinates 
   of a point in the @em reference (@c 2D) space
   @returns Returns the coordiantes of the given point
   in the @em physical space

   @fn Jacobian::map(const double, const double) const
   @param u The @c U coordinate 
   of a point in the @em reference (@c 2D) space
   @param v The @c V coordinate of the same point
   @returns Returns the coordiantes of the given point
   in the @em physical space
 */

//////////////////////
// Inline Functions //
//////////////////////

inline double Jacobian::det(void) const{
  return detDxDu;
}

#endif
