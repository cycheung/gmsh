#ifndef _GROUPOFELEMENT_H_
#define _GROUPOFELEMENT_H_

#include <string>
#include <vector>
#include <map>

#include "Mesh.h"
#include "MElement.h"

#include "Basis.h"
#include "Exception.h"

/**
   @class GroupOfElement
   @brief A Group of MElement%s

   This class is collection of @em discrete elements (MElement%s).
*/

class Mesh;

class GroupOfElement{
 private:
  const Mesh* mesh;

  std::vector<const MElement*>* element;
  std::vector<unsigned int>*    orientationStat;

 public:
   GroupOfElement(std::multimap<int, const MElement*>::iterator begin,
		  std::multimap<int, const MElement*>::iterator end,
		  const Mesh& mesh);

  ~GroupOfElement(void);

  unsigned int    getNumber(void)     const;
  const MElement& get(unsigned int i) const;

  const std::vector<const MElement*>&
    getAll(void) const;

  const Mesh& getMesh(void) const;

  void orientAllElements(const Basis& basis);
  const std::vector<unsigned int>& getOrientationStats(void) const;

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

   @fn GroupOfElement::orientAllElements
   @param basis A Basis

   Sort the Element of this GroupOfElement,
   with respect to the given Basis

   @see Basis::getOrientation()
   **

   @fn GroupOfElement::getOrientationStats
   @return A vector where the @c i-th entry is the number
   of element in GroupOfElement::getAll()
   with a Basis::getOrientation() equal to @c i

   @warning
   GroupOfElement::orientAllElement must be called
   before for this method to have a meaning@n

   If not, an Exception is thrown
   **

   @fn GroupOfElement::toString
   @return Returns a string discribing this Group
*/


//////////////////////
// Inline Functions //
//////////////////////

inline unsigned int GroupOfElement::getNumber(void) const{
  return element->size();
}

inline const MElement& GroupOfElement::get(unsigned int i) const{
  return *(*element)[i];
}

inline const std::vector<const MElement*>&
GroupOfElement::getAll(void) const{
  return *element;
}

inline const Mesh& GroupOfElement::getMesh(void) const{
  return *mesh;
}

inline const std::vector<unsigned int>&
GroupOfElement::getOrientationStats(void) const{
  if(orientationStat)
    return *orientationStat;

  else
    throw Exception("Orientation of GroupOfElement not computed");
}

#endif
