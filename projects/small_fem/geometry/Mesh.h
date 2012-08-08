#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <map>
#include <string>

#include "GroupOfElement.h"
#include "GModel.h"

/**
   @class Mesh
   @brief Represents a mesh
   
   This class represents a mesh.@n

   A Mesh is a collection of @em GroupOfElement%s 

   A Mesh is instantiated thanks to a 
   <a href="http://www.geuz.org/gmsh">gmsh</a>
   @c .msh file, wich discribes the mesh.@n
*/

class Mesh{
 private:
  GModel* model;

  unsigned int                        nEntity;
  std::vector<GroupOfElement*>*        group;
  std::multimap<int, GroupOfElement*>* physToGroup;

 public:
   Mesh(const std::string fileName);
  ~Mesh(void);

  unsigned int getGlobalId(const MVertex& vertex) const;
  unsigned int getGlobalId(const MEdge& edge) const;
  unsigned int getGlobalId(const MFace& face) const;
  unsigned int getGlobalId(const MElement& element) const;

  int                                 getNbGroup(void) const;
  GroupOfElement&                     getGroup(int i) const;
  const std::vector<GroupOfElement*>& getAllGroups(void) const;
  
  const std::vector<GroupOfElement*>  getFromPhysical(int physical) const;

  std::string toString(void) const;
};

/**
   @fn Mesh::Mesh
   Instanciate a new Mesh
   @param fileName The path to the @c .msh file discribing the Mesh
   
   @fn Mesh::~Mesh
   Deletes this Mesh

   @fn Mesh::getNbGroup
   @return Returns the number of group in this Mesh

   @fn Mesh::getGroup
   @param i A number between 0 and getNbGroup() - 1
   @return Returns the requested GroupOfElement

   @fn Mesh::getAllGroups
   @return Returns all the Mesh GroupOfElements
 
   @fn Mesh::toString
   @return Returns a description of this Mesh
*/


//////////////////////
// Inline Functions //
//////////////////////

inline unsigned int Mesh::getGlobalId(const MVertex& vertex) const{
  return vertex.getNum();
}

inline unsigned int Mesh::getGlobalId(const MElement& element) const{
  return element.getNum();
}

inline int Mesh::getNbGroup(void) const{
  return nEntity;
}

inline GroupOfElement& Mesh::getGroup(int i) const{
  return *((*group)[i]);
}

inline const std::vector<GroupOfElement*>&
Mesh::getAllGroups(void) const{
  return *group;
}

#endif
