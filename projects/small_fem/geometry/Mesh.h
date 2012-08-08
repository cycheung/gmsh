#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <map>
#include <string>

#include "GModel.h"
#include "GroupOfElement.h"
#include "GroupOfVertex.h"
#include "GroupOfEdge.h"

#include "MElement.h"
#include "MVertex.h"
#include "MEdge.h"
#include "MFace.h"

/**
   @class Mesh
   @brief Represents a mesh
   
   This class represents a mesh.@n

   A Mesh is a collection of @em GroupOfElement%s 

   A Mesh is instantiated thanks to a 
   <a href="http://www.geuz.org/gmsh">gmsh</a>
   @c .msh file, wich discribes the mesh.@n
*/

class GroupOfVertex;
class GroupOfEdge;

class MEdgeLessThanNum;
class MElementLessThanNum;

class Mesh{
 private:
  class MEdgeLessThanNum{
  public:
    bool operator()(const MEdge *e1, const MEdge *e2) const;
  };
  
  class MElementLessThanNum{
  public:
    bool operator()(const MElement *e1, const MElement *e2) const;
  };

 private:
  GModel* model;

  unsigned int                        nEntity;
  std::vector<GroupOfElement*>*        group;
  std::multimap<int, GroupOfElement*>* physToGroup;

  std::map<GroupOfElement*, 
           GroupOfVertex*, 
           GroupComparator>* elementToVertex;

  std::map<GroupOfElement*, 
           GroupOfEdge*, 
           GroupComparator>* elementToEdge;

  std::map<MVertex*, unsigned int, MVertexLessThanNum>*   vertex;
  std::map<MEdge*, unsigned int, MEdgeLessThanNum>*       edge;
  std::map<MElement*, unsigned int, MElementLessThanNum>* element;           

  std::map<unsigned int, MVertex*>*   idVertex;
  std::map<unsigned int, MEdge*>*       idEdge;
  std::map<unsigned int, MElement*>* idElement;           

  int nextEntityId;

 public:
   Mesh(const std::string fileName);
  ~Mesh(void);

  unsigned int getGlobalId(const MElement& element) const;
  unsigned int getGlobalId(const MVertex& vertex) const;
  unsigned int getGlobalId(const MEdge& edge) const;
  unsigned int getGlobalId(const MFace& face) const;

  MElement& getElement(unsigned int id) const;
  MVertex&  getVertex(unsigned int id) const;
  MEdge&    getEdge(unsigned int id) const;
  MFace&    getFace(unsigned int id) const;

  int                                 getNbGroup(void) const;
  GroupOfElement&                     getGroup(int i) const;
  const std::vector<GroupOfElement*>& getAllGroups(void) const;
  const std::vector<GroupOfElement*>  getFromPhysical(int physical) const;

  GroupOfVertex& getGroupOfVertex(GroupOfElement& goe);
  GroupOfEdge&   getGroupOfEdge(GroupOfElement& goe);
  
  std::string toString(void) const;

 private:
  static MEdge invert(MEdge& edge);

  void extractElement(GEntity& entity, int i);
  void extractVertex(GroupOfElement& goe);
  void extractEdge(GroupOfElement& goe);
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

inline bool Mesh::MEdgeLessThanNum::
operator()(const MEdge *e1, const MEdge *e2) const{
  if(e1->getMinVertex() < e2->getMinVertex()) return true;
  if(e1->getMinVertex() > e2->getMinVertex()) return false;
  if(e1->getMaxVertex() < e2->getMaxVertex()) return true;
  return false;
}

inline bool Mesh::MElementLessThanNum::
operator()(const MElement *e1, const MElement *e2) const{
  return e1->getNum() < e2->getNum();
}

#endif
