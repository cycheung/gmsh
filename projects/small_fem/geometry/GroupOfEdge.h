#ifndef _GROUPOFEDGE_H_
#define _GROUPOFEDGE_H_

#include <string>
#include <vector>

#include "Group.h"
#include "GroupOfElement.h"
#include "MEdge.h"
#include "MVertex.h"

/**
   @class GroupOfEdge
   @brief A Group of MEdge

   This class is collection of @em Edges (MEdge).@n
   This class is @em Group.
*/

class GroupOfElement;

class GroupOfEdge: public Group{
 private:
  const GroupOfElement* goe;

  unsigned int         nEdge;
  std::vector<MEdge*>*  edge;

 public:
  GroupOfEdge(const GroupOfElement& goe);
  virtual ~GroupOfEdge(void);

  virtual int getNumber(void) const;
  virtual int getType(void)   const;

  MEdge&                     get(int i) const;  
  const std::vector<MEdge*>& getAll(void) const;  

  virtual std::string toString(void) const;
};


//////////////////////
// Inline Functions //
//////////////////////

inline int GroupOfEdge::getNumber(void) const{
  return nEdge;
}

inline int GroupOfEdge::getType(void) const{
  return 3;
}

inline MEdge& GroupOfEdge::get(int i) const{
  return (*(*edge)[i]);
}

inline const std::vector<MEdge*>& 
GroupOfEdge::getAll(void) const{
  return *edge;
}

#endif
