#ifndef _GROUPOFEDGE_H_
#define _GROUPOFEDGE_H_

#include <string>
#include <vector>

#include "Group.h"
#include "GroupOfElement.h"
#include "MEdge.h"

/**
   @class GroupOfEdge
   @brief A Group of MEdge

   This class is collection of @em Edges (MEdge).@n
   This class is @em Group.
*/

class GroupOfEdge: public Group{
 public:
  GroupOfEdge(const GroupOfElement& goe);
  virtual ~GroupOfEdge(void);

  virtual int getNumber(void) const;
  virtual int getId(void)     const;
  virtual int getType(void)   const;

  MEdge&                     get(int i) const;  
  const std::vector<MEdge*>& getAll(void) const;  

  virtual std::string toString(void) const;
};

#endif
