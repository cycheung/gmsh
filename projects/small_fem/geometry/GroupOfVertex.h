#ifndef _GROUPOFVERTEX_H_
#define _GROUPOFVERTEX_H_

#include <string>
#include <vector>

#include "Mesh.h"
#include "Group.h"
#include "GroupOfElement.h"
#include "MVertex.h"

/**
   @class GroupOfVertex
   @brief A Group of MVertex

   This class is collection of @em Vertices (MVertex).@n
   This class is @em Group.
*/

class Mesh;
class GroupOfElement;

class GroupOfVertex: public Group{
 private:
  Mesh* mesh;

  static unsigned int nextId;
  unsigned int        id;

  unsigned int            nVertex;
  std::vector<MVertex*>*  vertex;

 public:
  GroupOfVertex(const GroupOfElement& goe, Mesh& mesh);
  virtual ~GroupOfVertex(void);

  virtual unsigned int getNumber(void) const;
  virtual unsigned int getId(void) const;
  virtual unsigned int getType(void)   const;

  MVertex&                     get(unsigned int i) const;  
  const std::vector<MVertex*>& getAll(void) const;  
  Mesh&                        getMesh(void) const;

  virtual std::string toString(void) const;
};


//////////////////////
// Inline Functions //
//////////////////////

inline unsigned int GroupOfVertex::getNumber(void) const{
  return nVertex;
}

inline unsigned int GroupOfVertex::getId(void) const{
  return id;
}

inline unsigned int GroupOfVertex::getType(void) const{
  return 2;
}

inline MVertex& GroupOfVertex::get(unsigned int i) const{
  return *((*vertex)[i]);
}

inline const std::vector<MVertex*>& 
GroupOfVertex::getAll(void) const{
  return *vertex;
}

inline Mesh& GroupOfVertex::getMesh(void) const{
  return *mesh;
}

#endif
