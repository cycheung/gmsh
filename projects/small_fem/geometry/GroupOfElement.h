#ifndef _GROUPOFELEMENT_H_
#define _GROUPOFELEMENT_H_

#include <string>
#include <vector>
#include <map>

#include "Mesh.h"
#include "MElement.h"


/**
   @class GroupOfElement
   @brief A Group of MElement%s

   This class is collection of @em discrete elements (MElement%s).
*/

class Mesh;

class GroupOfElement{
 private:
  const Mesh* mesh;
  
  unsigned int                  nElement;
  std::vector<const MElement*>*  element;

 public:
   GroupOfElement(std::multimap<int, const MElement*>::iterator begin, 
		  std::multimap<int, const MElement*>::iterator end,
		  const Mesh& mesh); 
 
  ~GroupOfElement(void);

  unsigned int                        getNumber(void)     const;
  const MElement&                     get(unsigned int i) const;  
  const std::vector<const MElement*>& getAll(void)        const;  
  
  const Mesh& getMesh(void) const;

  std::string toString(void) const;
};


/**
   @fn GroupOfElement::GroupOfElement
   @param begin An std::multimap @em Iterator
   @param end   An other std::mutltimap @em Iterator
   @param mesh  A Mesh

   Instantiates a new GroupOfElement, associated to the given Mesh
   **

   @fn GroupOfElement::~GroupOfElement
   Deletes this GroupOfElement
   **

   @fn GroupOfElement::getNumber
   @return Returns the number of elements in this GroupOfElement
   **

   @fn GroupOfElement::get
   @param i An interger ranging from 0 
   to GroupOfElement::getNumber() - 1
   @return Returns the ith element of the GroupOfElement
   **

   @fn GroupOfElement::getAll
   @return Returns all the elements of the GroupOfElement
   **

   @fn GroupOfElement::getMesh
   @return Returns the associated Mesh  
   **

   @fn GroupOfElement::toString
   @return Returns a string discribing this Group
*/


//////////////////////
// Inline Functions //
//////////////////////

inline unsigned int GroupOfElement::getNumber(void) const{
  return nElement;
}

inline const MElement& GroupOfElement::get(unsigned int i) const{
  return *((*element)[i]);
}

inline const std::vector<const MElement*>& 
GroupOfElement::getAll(void) const{
  return *element;
}

inline const Mesh& GroupOfElement::getMesh(void) const{
  return *mesh;
}

#endif
