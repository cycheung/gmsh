#ifndef _GROUPOFELEMENT_H_
#define _GROUPOFELEMENT_H_

#include <string>
#include <vector>
#include <map>

#include "Mesh.h"
#include "MElement.h"
#include "Basis.h"
#include "ElementType.h"

/**
   @class GroupOfElement
   @brief A Group of MElement%s

   This class is collection of discrete elements (MElement%s).
*/

class Mesh;

class GroupOfElement{
 private:
  class OrientationSort{
   private:
    const Basis* basis;

   public:
     OrientationSort(const Basis& basis);
    ~OrientationSort(void);

    bool operator()(const MElement* a, const MElement* b) const;
  };

 private:
  const Mesh* mesh;

  std::vector<const MElement*> element;
  std::vector<size_t>  orientationStat;

 public:
   GroupOfElement(std::multimap<int, const MElement*>::iterator begin,
                  std::multimap<int, const MElement*>::iterator end,
                  const Mesh& mesh);

  ~GroupOfElement(void);

  size_t                              getNumber(void) const;
  const MElement&                     get(size_t i)   const;
  const std::vector<const MElement*>& getAll(void)    const;
  const Mesh&                         getMesh(void)   const;

  void orientAllElements(const Basis& basis);
  const std::vector<size_t>& getOrientationStats(void) const;
  void unoriented(void);

  std::string toString(void) const;
};


/**
   @fn GroupOfElement::GroupOfElement
   @param begin An std::multimap Iterator
   @param end   An other std::mutltimap Iterator
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
   @param i An interger ranging from 0 to GroupOfElement::getNumber() - 1
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
   with respect to the given Basis orientations

   @see Basis::getOrientation()
   **

   @fn GroupOfElement::getOrientationStats
   @return A vector where the i-th entry is the number
   of element in GroupOfElement::getAll()
   with a Basis::getOrientation() equal to i

   GroupOfElement::orientAllElement must be called
   before for this method to have a meaning

   If not, an Exception is thrown
   **

   @fn GroupOfElement::unoriented

   The elements of this GroupOfElement are unoriented

   @see Basis::getOrientation()
   **

   @fn GroupOfElement::toString
   @return Returns a string discribing this Group
*/


//////////////////////
// Inline Functions //
//////////////////////

inline bool GroupOfElement::OrientationSort::operator()
(const MElement* a, const MElement* b) const{
  return
    basis->getReferenceSpace().getReferenceSpace(*a) <
    basis->getReferenceSpace().getReferenceSpace(*b);
}

inline const MElement& GroupOfElement::get(size_t i) const{
  return *element[i];
}

#endif
