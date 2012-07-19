#ifndef _GROUP_H_
#define _GROUP_H_

#include <vector>
#include <string>

/**
   @interface Group
   @brief A group of elements

   This is the interface allowing @em partial @em access to a group
   of elements.@n

   For a @em full access, a @em templated GroupTyped is requiered.
   @note Note that a GroupTyped is a Group.
 */

class Group{
 public:
  virtual ~Group(void);

  virtual int getNumber(void) const = 0;
  
  virtual int getId(void) const = 0; 

  virtual std::string toString(void) const = 0;
};

/**
   @fn Group::~Group
   Deletes this Group

   @fn Group::getNumberb
   @return Returns the number of elements in the Group

   @fn Group::getNbElements
   @return Returns the @c ID of this Group 

   @fn Group::toString
   @return Returns a string discribing this Group
*/

#endif
