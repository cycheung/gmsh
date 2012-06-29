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

   A Mesh is composed of @em physicals and of @em Group%s 
   @em of @em elements defined on those physicals.@n

   A Mesh is instantiated thanks to a 
   <a href="http://www.geuz.org/gmsh">gmsh</a>
   @c .msh file, wich discribes the mesh.@n
*/

class Mesh{
 private:
  GModel* model;
  
  std::map<int, Group*>* physToGroup;

 public:
   Mesh(const std::string fileName);
  ~Mesh(void);

  int    getNbPhysicals(void) const;
  Group& getGroup(int i) const;

  const std::vector<std::pair<int, Group*> >
    getAllGroups(void) const;
  
  std::string toString(void) const;
};

/**
   @fn Mesh::Mesh
   Instanciate a new Mesh
   @param fileName The path to the @c .msh file discribing the Mesh
   
   @fn Mesh::~Mesh
   Deletes this Mesh

   @fn Mesh::getNbPhysicals
   @return Returns the number of physical in this Mesh

   @fn Mesh::getGroup
   @param i A @em physical @c ID
   @return Returns the requested Group

   @fn Mesh::getAllGroups
   @return Returns all the pair <physical @c ID, Group>
 
   @fn Mesh::toString
   @return Returns a description of this Mesh
*/


//////////////////////
// Inline Functions //
//////////////////////

inline int Mesh::getNbPhysicals(void) const{
  return physToGroup->size();
}

inline Group& Mesh::getGroup(int i) const{
  return *(physToGroup->find(i)->second);
}

#endif
