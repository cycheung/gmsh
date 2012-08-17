#ifndef _GROUPOFELEMENT_H_
#define _GROUPOFELEMENT_H_

#include <string>
#include <vector>
#include <map>

#include "GroupTyped.h"
#include "Mesh.h"
#include "MElement.h"


/**
   @class GroupOfElement
   @brief A Group of MElement%s

   This class is collection (a Group) of @em discrete elements (MElement%s).@n
   
   In addition to the properties of a Group (GroupTyped<MElement>),@n
   a GroupOfElement can gives access the @em Mesh of its elements.
*/

class Mesh;

class GroupOfElement: public GroupTyped<MElement>{
 private:
  static unsigned int nextId;
  unsigned int        id;
  const Mesh*         mesh;
  
  unsigned int                  nElement;
  std::vector<const MElement*>*  element;

 public:
  GroupOfElement(std::multimap<int, const MElement*>::iterator begin, 
		 std::multimap<int, const MElement*>::iterator end,
		 const Mesh& mesh); 
 
  virtual ~GroupOfElement(void);

  virtual unsigned int getNumber(void) const;
  virtual unsigned int getId(void) const;
  virtual unsigned int getType(void)   const;

  virtual const MElement&                     get(unsigned int i) const;  
  virtual const std::vector<const MElement*>& getAll(void) const;  
  
  const Mesh& getMesh(void) const;

  virtual std::string toString(void) const;
};


/**
   @fn GroupOfElement::GroupOfElement
   @param begin An std::multimap @em Iterator
   @param end   An other std::mutltimap @em Iterator
   @param mesh  A Mesh

   Instantiates a new GroupOfElement, associated to the given Mesh@n

   The MElement%s of the Group are given by the two Iterators
   **

   @fn GroupOfElement::~GroupOfElement
   Deletes this GroupOfElement
   **

   @fn GroupOfElement::get
   @param i An interger ranging from 0 
   to GroupOfElement::getNumber() - 1
   @return Returns the ith element of the Group
   **

   @fn GroupOfElement::getAll
   @return Returns all the elements of the Group
   **

   @fn GroupOfElement::getMesh
   @return Returns the associated Mesh  
   **
*/


//////////////////////
// Inline Functions //
//////////////////////

inline unsigned int GroupOfElement::getNumber(void) const{
  return nElement;
}

inline unsigned int GroupOfElement::getId(void) const{
  return id;
}

inline unsigned int GroupOfElement::getType(void) const{
  return 1;
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
