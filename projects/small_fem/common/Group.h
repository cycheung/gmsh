#ifndef _GROUP_H_
#define _GROUP_H_

#include <vector>
#include <string>

/**
   @interface Group
   @brief A group of elements

   This is the interface allowing @em partial @em access to a group
   of elements.@n

   Every Group has a particular type of element.@n
   Each type is represented by a number:
   @li 0 for Group of @em Dof%s
   @li 1 for Group of @em MElement%s
   @li 2 for Group of @em MVertex%s
   @li 3 for Group of @em MEdge%s

   For a @em full access, use implemented classes (GroupOf*).@n

   @warning
   Every class implementing this interface @em must also implement:
   @li A get() method to access a @em particular element of the Group
   @li A getAll() method to access @em all the elements of the Group

   These methods are @em not in this interface because of @em return @em type
   @em issues.

   @todo
   Find a solution for return type issuses, template and doxygen !!
 */

/**
   @class GroupComparator
   @brief A class to compare two Group%s

   A class to compare two Group%s of the @em same type.
*/


class Group{
 public:
  virtual ~Group(void);

  virtual unsigned int getNumber(void) const = 0;
  virtual unsigned int getId(void) const = 0;
  virtual unsigned int getType(void) const = 0;

  virtual std::string toString(void) const = 0;
};

class GroupComparator{
 public:
  bool operator()(const Group* a, 
		  const Group* b) const;
};

/**
   @fn Group::~Group
   Deletes this Group

   @fn Group::getNumber
   @return Returns the number of elements in the Group

   @fn Group::getId
   @return Returns the (unique) @c ID of this Group
   @note 
   An @c ID is unique @em for @em a @em given @em type of Group 

   @fn Group::getType
   @return Returns the type of the elements 

   @fn Group::toString
   @return Returns a string discribing this Group

   @fn bool GroupComparator::operator(Group* a, const Group* b) const
   @param a A Group
   @param b Another Group
   @return operator() is:
   @li @c true, if a is @em smaller than b  
   @li @c false, otherwise
*/


//////////////////////
// Inline Functions //
//////////////////////

inline bool GroupComparator::operator()(const Group* a, 
					const Group* b) const{
  return a->getId() < b->getId();
}

#endif
