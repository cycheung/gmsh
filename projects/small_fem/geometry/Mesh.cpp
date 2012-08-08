#include "Mesh.h"
#include "Exception.h"

#include <list>
#include <sstream>

using namespace std;

Mesh::Mesh(const std::string fileName){ 
  // New Mode //
  model = new GModel("SmallFEM");

  // Read Mesh //
  if(!model->readMSH(fileName))
    throw Exception("Can't open file: %s", fileName.c_str());

  // Get Entity //
  vector<GEntity*> entity;
  model->getEntities(entity);
  nEntity = entity.size();

  // Alloc Memory //
  group       = new vector<GroupOfElement*>(nEntity);
  physToGroup = new multimap<int, GroupOfElement*>;
  
  elementToVertex = new map<GroupOfElement*,
			    GroupOfVertex*,
			    GroupComparator>;

  elementToEdge   = new map<GroupOfElement*,
			    GroupOfEdge*,
			    GroupComparator>;

  vertex  = new map<MVertex*, unsigned int, MVertexLessThanNum>;
  edge    = new map<MEdge*, unsigned int, MEdgeLessThanNum>;
  element = new map<MElement*, unsigned int, MElementLessThanNum>;

  idVertex  = new map<unsigned int, MVertex*>;
  idEdge    = new map<unsigned int, MEdge*>;
  idElement = new map<unsigned int, MElement*>;

  nextEntityId = 0;

  // Extract Element
  for(unsigned int i = 0; i < nEntity; i++)
    extractElement(*entity[i], i);

  // Extract Vertex
  for(unsigned int i = 0; i < nEntity; i++)
    extractVertex(*(*group)[i]);

  // Extract Edge
  for(unsigned int i = 0; i < nEntity; i++)
    extractEdge(*(*group)[i]);
}

Mesh::~Mesh(void){
  // Iterator //
  map<GroupOfElement*,
      GroupOfVertex*,
      GroupComparator>::iterator itEtV;

  map<GroupOfElement*,
      GroupOfEdge*,
      GroupComparator>::iterator itEtE;

  // End Iterator //
  const map<GroupOfElement*,
	    GroupOfVertex*,
	    GroupComparator>::iterator endEtV = 
    elementToVertex->end();

  const map<GroupOfElement*,
	    GroupOfEdge*,
	    GroupComparator>::iterator endEtE = 
    elementToEdge->end();

  // Delete Element To Vertex //
  for(itEtV = elementToVertex->begin(); 
      itEtV != endEtV; itEtV++)
    delete itEtV->second;

  delete elementToVertex;

  // Delete Element To Edge //
  for(itEtE = elementToEdge->begin(); 
      itEtE != endEtE; itEtE++)
    delete itEtE->second;

  delete elementToEdge;

  // Delete Vertex //
  delete vertex;
  delete idVertex;

  // Delete Edge //
  delete edge;
  delete idEdge;

  // Delete Element //
  delete element;
  delete idElement;

  // Delete Model //
  delete model;  

  // Delete Phys To Group //
  delete physToGroup;

  // Delete Group //
  for(unsigned int i = 0; i < nEntity; i++)
    delete (*group)[i];
  delete group;
}

const vector<GroupOfElement*> Mesh::getFromPhysical(int physical) const{
  pair<multimap<int, GroupOfElement*>::iterator, 
       multimap<int, GroupOfElement*>::iterator> startStop = 
    physToGroup->equal_range(physical);
  
  multimap<int, GroupOfElement*>::iterator it;
  list<GroupOfElement*> lst;
  
  for(it = startStop.first; it != startStop.second; it++)
    lst.push_back((*it).second);

  return vector<GroupOfElement*>(lst.begin(), lst.end());
}

unsigned int Mesh::getGlobalId(const MElement& element) const{
  MElement& eelement = const_cast<MElement&>(element);

  map<MElement*, 
      unsigned int, 
      MElementLessThanNum>::iterator it = 
    this->element->find(&eelement);

  if(it == this->element->end())
    throw 
      Exception("Element not found");

  return it->second;   
}

unsigned int Mesh::getGlobalId(const MVertex& vertex) const{
  MVertex& vvertex = const_cast<MVertex&>(vertex);

  map<MVertex*, 
      unsigned int, 
      MVertexLessThanNum>::iterator it = 
    this->vertex->find(&vvertex);

  if(it == this->vertex->end())
    throw 
      Exception("Vertex not found");

  return it->second; 
}

unsigned int Mesh::getGlobalId(const MEdge& edge) const{
  MEdge& eedge = const_cast<MEdge&>(edge);

  // WARNING:
  // Here, we use a Edge Comparator,
  // such that Edge Orientation
  // do *NOT* matter !!

  // Look for Edge //
  map<MEdge*, 
      unsigned int, 
      MEdgeLessThanNum>::iterator it = 
    this->edge->find(&eedge);

  if(it == this->edge->end()){
    throw 
      Exception("Edge not found");
  }

  return it->second; 
}

MElement& Mesh::getElement(unsigned int id) const{
  map<unsigned int, MElement*>::iterator it = 
    idElement->find(id);

  if(it == idElement->end())
    throw 
      Exception("No Element with Global Id %d found", id);

  return *(it->second);   
}

MVertex& Mesh::getVertex(unsigned int id) const{
  map<unsigned int, MVertex*>::iterator it = 
    idVertex->find(id);

  if(it == idVertex->end())
    throw 
      Exception("No Vertex with Global Id %d found", id);

  return *(it->second);   
}

MEdge& Mesh::getEdge(unsigned int id) const{
  map<unsigned int, MEdge*>::iterator it = 
    idEdge->find(id);

  if(it == idEdge->end())
    throw 
      Exception("No Edge with Global Id %d found", id);

  return *(it->second);   
}

GroupOfVertex& Mesh::getGroupOfVertex(GroupOfElement& goe){
  // Iterator //
  map<GroupOfElement*, 
      GroupOfVertex*, 
      GroupComparator>::iterator it = elementToVertex->find(&goe);
  
  map<GroupOfElement*, 
      GroupOfVertex*, 
      GroupComparator>::iterator end = elementToVertex->end();

  // Check If we have the requested Group
  if(it != end)
    return *(it->second);
  
  // Else: Exception !
  else
    throw Exception("GroupOfVertex not found");
}

GroupOfEdge& Mesh::getGroupOfEdge(GroupOfElement& goe){
  // Iterator //
  map<GroupOfElement*, 
      GroupOfEdge*, 
      GroupComparator>::iterator it = elementToEdge->find(&goe);
  
  map<GroupOfElement*, 
      GroupOfEdge*, 
      GroupComparator>::iterator end = elementToEdge->end();

  // Check If we have the requested Group
  if(it != end)
    return *(it->second);
  
  // Else: Exception !
  else
    throw Exception("GroupOfVertex not found");
}

void Mesh::extractElement(GEntity& entity, int i){
  // Build & Add GroupOfElement
  GroupOfElement* goe = new GroupOfElement(entity, *this); 
  (*group)[i]         = goe;

  // Add physical value
  vector<int> physical = entity.getPhysicalEntities();
  int nPhysical        = physical.size();

  for(int j = 0; j < nPhysical; j++)
    physToGroup->insert(pair<int, GroupOfElement*>(physical[j], goe));

  // Number Element
  const std::vector<MElement*>& myElement = goe->getAll();
  unsigned int                          N = goe->getNumber();
  
  for(unsigned int j = 0; j < N; j++){
    element->insert
      (pair<MElement*, unsigned int>
       (myElement[j], nextEntityId));
    
    idElement->insert
      (pair<unsigned int, MElement*>
       (nextEntityId, myElement[j]));
    
    nextEntityId++;
  }  
}

void Mesh::extractVertex(GroupOfElement& goe){
  // Build GroupOfVertex
  GroupOfVertex* gov = new GroupOfVertex(goe, *this);
  
  // Insert in lookup table
  elementToVertex->insert
    (pair<GroupOfElement*, GroupOfVertex*>(&goe, gov));

  // Number Vertices
  const std::vector<MVertex*>& myVertex = gov->getAll();
  unsigned int                        N = gov->getNumber();
  
  for(unsigned int i = 0; i < N; i++){
    vertex->insert
      (pair<MVertex*, unsigned int>
       (myVertex[i], nextEntityId));
    
    idVertex->insert
      (pair<unsigned int, MVertex*>
       (nextEntityId, myVertex[i]));
    
    nextEntityId++;
  }
}

void Mesh::extractEdge(GroupOfElement& goe){
  // Build GroupOfEdge
  GroupOfEdge* goEd = new GroupOfEdge(goe, *this);
  
  // Insert in lookup table
  elementToEdge->insert
    (pair<GroupOfElement*, GroupOfEdge*>(&goe, goEd));

  // Number Edges
  const std::vector<MEdge*>& myEdge = goEd->getAll();
  unsigned int                    N = goEd->getNumber();
  
  for(unsigned int i = 0; i < N; i++){
    edge->insert
      (pair<MEdge*, unsigned int>
       (myEdge[i], nextEntityId));
    
    idEdge->insert
      (pair<unsigned int, MEdge*>
       (nextEntityId, myEdge[i]));
    
    nextEntityId++;
  }
}

string Mesh::toString(void) const{
  stringstream stream;

  stream << "*********************************************"    
	 << endl
	 << "*                    Mesh                   *"    
	 << endl
	 << "*********************************************"
	 << endl << endl;
  

  for(unsigned int i = 0; i < nEntity; i++)
    stream << (*group)[i]->toString() << endl;

  return stream.str();
}

MEdge Mesh::invert(MEdge& edge){
  MVertex* begin = edge.getVertex(0);
  MVertex* end   = edge.getVertex(1);

  return MEdge(end, begin);
}
