#ifndef _GROUP_H_
#define _GROUP_H_

#include <string>

/**
   @interface Group
   @brief A group of elements (with partial access)

   This is the interface allowing @em partial @em access to a group
   of elements.@n

   Every Group has a particular type of element.@n
   Each type is represented by a number:
   @li 0 for Group of @em Dof%s
   @li 1 for Group of @em MElement%s
   @li 2 for Group of @em MVertex%s
   @li 3 for Group of @em MEdge%s

   A Group (for a  @em given @em type of element) 
   got also a @em unique @c ID.

   For a @em full access, use GroupTyped (with required template).@n


   @note
   Note that a GroupTyped is a Group
 */

class Group{
 public:
  virtual ~Group(void);

  virtual unsigned int getNumber(void) const = 0;
  virtual unsigned int getId(void) const = 0;
  virtual unsigned int getType(void) const = 0;

  virtual std::string toString(void) const = 0;
};


/**
   @fn Group::~Group
   Deletes this Group
   **

   @fn Group::getNumber
   @return Returns the number of elements in the Group
   **

   @fn Group::getId
   @return Returns the (unique) @c ID of this Group
   @note 
   An @c ID is unique @em for @em a @em given @em type of Group 
   **

   @fn Group::getType
   @return Returns the type of the elements 
   **

   @fn Group::toString
   @return Returns a string discribing this Group
   **
*/


#endif
