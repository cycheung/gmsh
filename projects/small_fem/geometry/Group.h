#ifndef _GROUP_H_
#define _GROUP_H_

#include <vector>
#include <string>
#include "MElement.h"

/**
   @interface Group
   @brief A group of elements

   This is the interface allowing @em access to a group
   of @em mesh (discrete) elements.
 */

class Group{
 public:
  virtual ~Group(void);

  virtual int getNumElements(void) const = 0;
  
  virtual MElement& getElement(int i) const = 0;
  
  virtual const std::vector<MElement*> getAllElements(void) const = 0;

  virtual std::string toString(void) const;
};

/**
   @fn Group::~Group
   Deletes this Group

   @fn Group::getNumElements
   @return Returns the number of elements in the Group
 
   @fn Group::getElement
   @param i An index ranging from 0 to getNumElements() - 1
   @return Returns the ith element of the Group

   @fn Group::getAllElements
   @return Returns all the elements of the Group

   @fn Group::toString
   @return Returns a string discribing this Group
*/

#endif
