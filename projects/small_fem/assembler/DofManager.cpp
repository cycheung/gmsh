#include <sstream>
#include "GroupOfElement.h"
#include "MVertex.h"
#include "MEdge.h"
#include "MFace.h"
#include "DofManager.h"



#include <iostream>

using namespace std;

DofManager::DofManager(const FunctionSpace& fs){
  // Remember FunctionSpace //
  this->fs = &fs;

  // Get Support from FunctionSpace //
  const GroupOfElement& support    = fs.getSupport();
  int nElement                     = support.getNumber();
  const vector<MElement*>& element = support.getAll();
  
  nTotVertex = support.getNVertex();

  // Init Struct //
  dof      = new set<Dof*, DofComparator>;         
  group    = new vector<GroupOfDof*>(nElement);
  globalId = new map<Dof*, int, DofComparator>;

  // Create Dofs & Numbering//
  nextId = 0;

  for(int i = 0; i < nElement; i++)
    add(*(element[i]), i);
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
  const int nDofVertex = nFVertex * nVertex; 
  const int nDofEdge   = nFEdge   * nEdge;
  const int nDofFace   = nFFace   * nFace;
  const int nDofCell   = nFCell;

  const int nDof = 
    nDofVertex + nDofEdge + nDofFace + nDofCell;
  
  (*group)[groupId] = new GroupOfDof(nDof, element);
  
  // Add Vertex Based Dof //
  for(int i = 0; i < nVertex; i++){
    // Get Id of Vertex
    const int id = vertex[i]->getNum();

    // Insert new Dof
    for(int j = 0; j < nFVertex; j++){
      cout << "Inserting Vertex (" << id << "): ";

      Dof* tmp = new Dof(id, j);
      insertDof(tmp, (*group)[groupId]);
    }
  }

  // Add Edge Based Dof //
  for(int i = 0; i < nEdge; i++){
    // Get Id of Edge 
    MVertex* vEdge0 = edge[i].getSortedVertex(0);
    MVertex* vEdge1 = edge[i].getSortedVertex(1);

    const int id = 
      vEdge0->getNum() + 
      vEdge1->getNum() * nTotVertex;

    // Insert new Dof
    for(int j = 0; j < nFEdge; j++){
      cout << "Inserting Edge" 
	   << "(" << vEdge0->getNum() << ", " << vEdge1->getNum() << ") "
	   << "(" << id << "): ";
      
      Dof* tmp = new Dof(id, j);
      insertDof(tmp, (*group)[groupId]);
    }
  }  
}

void DofManager::insertDof(Dof* d, GroupOfDof* god){
  // Try to insert Dof //
  pair<set<Dof*, DofComparator>::iterator, bool> p;
  p = dof->insert(d);
 
  // If insertion is OK (Dof 'd' didn't exist) //
  //   --> Add new Dof
  if(p.second){
    cout << "Yes -- ID: " << nextId << endl;

    globalId->insert(pair<Dof*, int>(d, nextId));
    
    god->add(d);
    
    nextId += 1;
  }
  
  // If insertion failed (Dof 'd' already exists) //
  //   --> delete 'd' and add existing Dof
  else{
    cout << "No" << endl;

    delete d; 
    god->add(*(p.first));
  }
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
