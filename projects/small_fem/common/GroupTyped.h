#ifndef _GROUPTYPED_H_
#define _GROUPTYPED_H_

#include <vector>
#include <string>
#include "Group.h"

/**
   @interface GroupTyped
   @brief A Group of typed elements

   This is the interface allowing @em full @em access to a group
   of elements.@n

   This class is templated according to the @em type of the stored 
   elements.
 */

template <class T>
class GroupTyped: public Group{
 public:
  virtual ~GroupTyped(void);

  virtual T& get(int i) const = 0;
  
  virtual const std::vector<T*>& getAll(void) const = 0;
};

/**
   @fn GroupTyped::~GroupTyped
   Deletes this GroupTyped
 
   @fn GroupTyped::get
   @param i An index ranging from 0 to getNumber() - 1
   @return Returns the ith element of the GroupTyped

   @fn GroupTyped::getAll
   @return Returns all the elements of the GroupTyped
*/


/////////////////////////
// Templated Functions //
/////////////////////////

template <class T>
GroupTyped<T>::~GroupTyped(void){
}

#endif
