#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <map>
#include <string>

#include "Group.h"
#include "GModel.h"

/**
   @class Mesh
   @brief Represents a mesh
   
   This class represents a mesh.@n

   A Mesh is a collection of @em Group%s 

   A Mesh is instantiated thanks to a 
   <a href="http://www.geuz.org/gmsh">gmsh</a>
   @c .msh file, wich discribes the mesh.@n
*/

class Mesh{
 private:
  GModel* model;

  unsigned int               nEntity;
  std::vector<Group*>*        group;
  std::multimap<int, Group*>* physToGroup;

 public:
   Mesh(const std::string fileName);
  ~Mesh(void);

  int                        getNbGroup(void) const;
  Group&                     getGroup(int i) const;
  const std::vector<Group*>& getAllGroups(void) const;
  
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
   @return Returns the requested Group

   @fn Mesh::getAllGroups
   @return Returns all the Mesh Groups
 
   @fn Mesh::toString
   @return Returns a description of this Mesh
*/


//////////////////////
// Inline Functions //
//////////////////////

inline int Mesh::getNbGroup(void) const{
  return nEntity;
}

inline Group& Mesh::getGroup(int i) const{
  return *((*group)[i]);
}

inline const std::vector<Group*>&
Mesh::getAllGroups(void) const{
  return *group;
}

#endif
