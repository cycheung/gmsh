#ifndef _GEODOF_H_
#define _GEODOF_H_

#include <vector>
#include "Dof.h"
#include "Mapper.h"
#include "MElement.h"
#include "fullMatrix.h"

/**
   @class GeoDof
   @brief Handel a group of Dof%s with @em geometrical meaning

   This class handles a group of Dof%s with a @em geometrical meaning 
   (@e e.g: Dof%s that belongs to the same (finite) element).@n

   This class allows geometrical computations (Reference to Physical
   space mapping, orientations, etc).

   @todo
   Use 'const MElement' --> need to change MElement::*Jacobian with const qualifire
*/


class DofManager;

class GeoDof{
 private:
  MElement* element;

  int nDof;
  std::vector<Dof*>* dof;
  const std::vector<int>* direction;
  
  int nextDof;
  
  friend class DofManager;

 public:
  int getId(void) const;
  int dofNumber(void) const;
  const std::vector<Dof*>& getAllDofs(void) const;

  int getOrientation(const int dofId) const;

  double        getJacobian(double u, double v, double w)      const;
  fullVector<double> map   (const fullVector<double>& UVW)     const; 
  fullVector<double> grad  (const fullVector<double>& gradUVW) const;
  fullVector<double> invMap(const fullVector<double>& XYZ)     const;

 private:
   GeoDof(int numberOfDof, const MElement& geoElement);
  ~GeoDof(void);

  void add(Dof* dof);
  void orientation(const std::vector<int>& orientation);
};


/**
   @fn int GeoDof::getId(void) const
   @return Returns the @c ID of this GeoDof

   @fn int GeoDof::dofNumber(void) const
   @return Returns the number of Dof in this GeoDof

   @fn const std::vector<Dof*>& GeoDof::getAllDofs(void) const
   @return Returns all the Dof%s in this GeoDof

   @fn const Jacobian& GeoDof::getJacobian(void) const;
   @return Returns the Jacobian associated to this GeoDof
   
   @fn int GeoDof::getOrientation(const int dofId) const;
   @param dofId The @em local @c ID of a Dof in the GeoDof
   @return Returns the orientation of a Dof (indentified by its @em local @c ID)
*/

//////////////////////
// Inline Functions //
//////////////////////

inline void GeoDof::orientation(const std::vector<int>& orientation){
  direction = &orientation;
}

inline int GeoDof::getId(void) const{
  return element->getNum();
}

inline int GeoDof::dofNumber(void) const{
  return nDof;
}

inline const std::vector<Dof*>& GeoDof::getAllDofs(void) const{
  return *dof;
}

inline int GeoDof::getOrientation(const int dofId) const{
  return (*direction)[dofId];
}


inline double GeoDof::getJacobian(double u, double v, double w) const{
  return Mapper::det(u, v, w, *(element));
}

inline fullVector<double> GeoDof::map(const fullVector<double>& UVW) const{
  return Mapper::map(UVW, *(element));
}

inline fullVector<double> GeoDof::grad(const fullVector<double>& gradUVW) const{
  return Mapper::grad(gradUVW, *(element));
}

inline fullVector<double> GeoDof::invMap(const fullVector<double>& XYZ) const{
  return Mapper::invMap(XYZ, *(element));
}

#endif
