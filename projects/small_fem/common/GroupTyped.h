#ifndef _GROUPTYPED_H_
#define _GROUPTYPED_H_

#include <vector>
#include "Group.h"

/**
   @interface GroupTyped
   @brief A Group of elements (with full access)

   This is the interface allowing @em full @em access 
   to a member of Group.@n

   This mean we can access all the element of the the Group.@n

   @note
   Note that a GroupTyped is a Group
 */

template<class T>
class GroupTyped: public Group{
 public:
  virtual ~GroupTyped(void);

  virtual const T& get(unsigned int i) const = 0; 
  virtual const std::vector<const T*>& getAll(void) const = 0;
};

/**
   @fn GroupTyped::~GroupTyped
   Deletes this GroupTyped
   **

   @fn GroupTyped::get
   @param i An interger ranging from 0 
   to getNumber() - 1
   @return Returns the ith element of the Group
   **

   @fn GroupTyped::getAll
   @return Returns all the elements of the Group
   **
*/


/////////////////////////
// Templated Functions //
/////////////////////////

template<class T>
GroupTyped<T>::~GroupTyped(void){
}

#endif
