#ifndef _GROUPOFVERTEX_H_
#define _GROUPOFVERTEX_H_

#include <string>
#include <vector>

#include "Group.h"
#include "GroupOfElement.h"
#include "MVertex.h"

/**
   @class GroupOfVertex
   @brief A Group of MVertex

   This class is collection of @em Vertices (MVertex).@n
   This class is @em Group.
*/

class GroupOfElement;

class GroupOfVertex: public Group{
 private:
  unsigned int            nVertex;
  std::vector<MVertex*>*  vertex;

 public:
  GroupOfVertex(const GroupOfElement& goe);
  virtual ~GroupOfVertex(void);

  virtual int getNumber(void) const;
  virtual int getType(void)   const;

  MVertex&                     get(int i) const;  
  const std::vector<MVertex*>& getAll(void) const;  

  virtual std::string toString(void) const;
};


//////////////////////
// Inline Functions //
//////////////////////

inline int GroupOfVertex::getNumber(void) const{
  return nVertex;
}

inline int GroupOfVertex::getType(void) const{
  return 2;
}

inline MVertex& GroupOfVertex::get(int i) const{
  return *((*vertex)[i]);
}

inline const std::vector<MVertex*>& 
GroupOfVertex::getAll(void) const{
  return *vertex;
}

#endif
