#ifndef _GROUPOFDOF_H_
#define _GROUPOFDOF_H_

#include <vector>

#include "Dof.h"
#include "MElement.h"

/**
   @class GroupOfDof
   @brief Handels a Group of Dof%s with geometrical meaning

   This class handles a collection of Dof%s that are
   associated to the same mesh Element.

   It also gives acces to the underlying mesh Element.
*/

class GroupOfDof{
 private:
  const MElement*  element;
  std::vector<Dof> dof;

 public:
  GroupOfDof(const MElement& element,
             const std::vector<Dof>& dof);

  ~GroupOfDof(void);

  size_t                  size(void)       const;
  const std::vector<Dof>& getDof(void)     const;
  const MElement&         getElement(void) const;

  std::string toString(void) const;
};


/**
   @fn GroupOfDof::GroupOfDof
   @param element A mesh Element (MElement)
   @param dof A vector of Dof

   Instanciates a new GroupOfDof related to the given
   Element and that can contains the given Dof%s
   **

   @fn GroupOfDof::~GroupOfDof

   Deletes this GroupOfDof
   **

   @fn GroupOfDof::size
   @return Returns the number of Dof%s in this GroupOfDof
   **

   @fn GroupOfDof::getDof
   @return Returns all the Dof%s
   **

   @fn GroupOfDof::getElement
   @return Returns the underlying mesh Element
   **

   @fn GroupOfDof::toString
   @return Returns a string discribing this GroupOfDof
*/

//////////////////////
// Inline Functions //
//////////////////////

inline size_t GroupOfDof::size(void) const{
  return dof.size();
}

inline const std::vector<Dof>& GroupOfDof::getDof(void) const{
  return dof;
}

inline const MElement& GroupOfDof::getElement(void) const{
  return *element;
}

#endif
