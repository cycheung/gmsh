#include <sstream>
#include "GroupOfElement.h"
#include "MVertex.h"
#include "MEdge.h"
#include "MFace.h"
#include "DofManager.h"

using namespace std;

DofManager::DofManager(const FunctionSpace& fs){
  // Remember FunctionSpace //
  this->fs = &fs;

  // Get Support from FunctionSpace //
  const GroupOfElement& support    = fs.getSupport();
  int nElement                     = support.getNumber();
  const vector<MElement*>& element = support.getAll();

  // Init Struct //
  dof      = new set<Dof*>;         
  group    = new vector<GroupOfDof*>(nElement);
  globalId = new map<Dof*, int, DofComparator>;

  // Create Dofs & Numbering//
  nextId = 0;

  for(int i = 0; i < nElement; i++)
    add(*(element[i]), i);
}


DofManager::DofManager(void){
}

DofManager::~DofManager(void){
  int nElement = group->size();

  for(int i = 0; i < nElement; i++)
    delete (*group)[i];
  delete group;

  set<Dof*>::iterator it;
  set<Dof*>::iterator end = dof->end();

  for(it = dof->begin(); it != end; it++)
    delete *it;
  delete dof;

  delete globalId;
}

void DofManager::add(MElement& element, int groupId){  
  // Get Element Data //
  const int nVertex = element.getNumVertices();
  const int nEdge   = element.getNumEdges();
  const int nFace   = element.getNumFaces(); 

  vector<MVertex*> vertex(nVertex);
  vector<MEdge> edge(nEdge);
  vector<MFace> face(nFace);

  for(int i = 0; i < nVertex; i++)
    vertex[i] = element.getVertex(i);

  for(int i = 0; i < nEdge; i++)    
    edge[i] = element.getEdge(i);
  
  for(int i = 0; i < nFace; i++)
    face[i] = element.getFace(i);
  
  // Get FunctionSpace Data for this Element //
  const int nFVertex = fs->getNFunctionPerVertex(element);
  const int nFEdge   = fs->getNFunctionPerEdge(element);
  const int nFFace   = fs->getNFunctionPerFace(element);
  const int nFCell   = fs->getNFunctionPerCell(element);

  // Create GroupOfDof //
  const int nDof = 
    nFVertex * nVertex +
    nFEdge   * nEdge   +
    nFFace   * nFace   +
    nFCell;
  
  (*group)[groupId] = new GroupOfDof(nDof, element);
  /*
  for(int i = 0; i < nEntity; i++){
    pair<set<Dof*, DofComparator>::iterator, bool> p;
    Dof* tmp = new Dof(entity[i]->getNum(), type);

    p = dofLookup->insert(tmp);
 
    if(p.second){
      globalId->insert(pair<Dof*, int>(tmp, nextId));
      
      (*group)[groupId]->add(tmp);

      nextId += 1;
    }

    else{
      delete tmp; // Dof already exists
      (*group)[groupId]->add(*(p.first)); // Add real Dof
    }
  }
  */
}

string DofManager::toString(void) const{
  stringstream s;
  map<Dof*, int, DofComparator>::iterator i   = globalId->begin();
  map<Dof*, int, DofComparator>::iterator end = globalId->end();
  
  for(; i != end; i++)
    s << "("  << (*i).first->toString() 

      << ": " << (*i).second 
      << ")"  << endl;

  return s.str();
}
