#ifndef _GROUPOFDOF_H_
#define _GROUPOFDOF_H_

#include <vector>

#include "Dof.h"
#include "MElement.h"

/**
   @class GroupOfDof
   @brief Handels a Group of Dof%s with @em geometrical meaning

   This class handles a collection of Dof%s with a @em geometrical meaning
   (@e e.g: Dof%s that belongs to the same (finite) element).@n

   It also gives acces to the underlying Geometrical Element.
*/

class GroupOfDof{
 private:
  const MElement*          element;
  std::vector<const Dof*>* dof;

 public:
  GroupOfDof(const MElement& geoElement,
             const std::vector<const Dof*>& dof);

  ~GroupOfDof(void);

  size_t                         size(void)          const;
  const std::vector<const Dof*>& getDof(void)        const;
  const MElement&                getGeoElement(void) const;

  std::string toString(void) const;
};


/**
   @fn GroupOfDof::GroupOfDof
   @param geoElement A geomtrical Element (MElement)
   @param dof A vector of Dof

   Instanciates a new GroupOfDof related to the given
   Element and that can contains the given Dof%s
   **

   @fn GroupOfDof::~GroupOfDof

   Deletes this GroupOfDof
   **

   @fn GroupOfDof::size
   @return Returns the number of elements in this GroupOfDof
   **

   @fn GroupOfDof::getDof
   @return Returns all the Dofs
   **

   @fn GroupOfDof::getGeoElement
   @return Returns the underlying Geometrical Element
   **

   @fn GroupOfDof::toString
   @return Returns a string discribing this GroupOfDof
*/

//////////////////////
// Inline Functions //
//////////////////////

inline size_t GroupOfDof::size(void) const{
  return dof->size();
}

inline const std::vector<const Dof*>& GroupOfDof::getDof(void) const{
  return *dof;
}

inline const MElement& GroupOfDof::getGeoElement(void) const{
  return *element;
}

#endif
