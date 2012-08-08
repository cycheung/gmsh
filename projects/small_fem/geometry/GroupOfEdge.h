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
  static unsigned int nextId;
  unsigned int            id;

  const GroupOfElement* goe;

  unsigned int         nEdge;
  std::vector<MEdge*>*  edge;

 public:
  GroupOfEdge(const GroupOfElement& goe);
  virtual ~GroupOfEdge(void);

  virtual unsigned int getNumber(void) const;
  virtual unsigned int getId(void) const;
  virtual unsigned int getType(void)   const;

  MEdge&                     get(unsigned int i) const;  
  const std::vector<MEdge*>& getAll(void) const;  

  virtual std::string toString(void) const;
};


//////////////////////
// Inline Functions //
//////////////////////

inline unsigned int GroupOfEdge::getNumber(void) const{
  return nEdge;
}

inline unsigned int GroupOfEdge::getId(void) const{
  return id;
}

inline unsigned int GroupOfEdge::getType(void) const{
  return 3;
}

inline MEdge& GroupOfEdge::get(unsigned int i) const{
  return (*(*edge)[i]);
}

inline const std::vector<MEdge*>& 
GroupOfEdge::getAll(void) const{
  return *edge;
}

#endif
