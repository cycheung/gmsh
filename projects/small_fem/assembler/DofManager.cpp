#include <sstream>
#include "GroupOfElement.h"
#include "MVertex.h"
#include "MEdge.h"
#include "MFace.h"
#include "Exception.h"
#include "DofManager.h"

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
  globalId = new map<const Dof*, int, DofComparator>;
  fixedDof = new map<const Dof*, double, DofComparator>;

  // Create Dofs & Numbering//
  nextId = 0;

  // Loop over Element
  for(int i = 0; i < nElement; i++){
    // Get Dof for this Element
    int nDof;
    Dof** myDof = dofFromElement(*(element[i]) ,&nDof);
    
    // Create new GroupOfDof
    (*group)[i] = new GroupOfDof(nDof, *(element[i]));

    // Add Dof
    for(int j = 0; j < nDof; j++)
      insertDof(myDof[j], (*group)[i]);

    // Delete Dof 'Container'
    delete[] myDof;
  }
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
  delete fixedDof;
}

Dof** DofManager::dofFromElement(MElement& element, int* nDof){  
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

  // Create Dof //
  const int nDofVertex = nFVertex * nVertex; 
  const int nDofEdge   = nFEdge   * nEdge;
  const int nDofFace   = nFFace   * nFace;
  const int nDofCell   = nFCell;

  *nDof = 
    nDofVertex + nDofEdge + nDofFace + nDofCell;

  Dof** myDof = new Dof*[*nDof];
  
  int it = 0;

  // Add Vertex Based Dof //
  for(int i = 0; i < nVertex; i++){
    // Get Id of Vertex
    const int id = vertex[i]->getNum();

    // New Dof
    for(int j = 0; j < nFVertex; j++){
      myDof[it] = new Dof(id, j);
      it++;
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
      myDof[it] = new Dof(id, j);
      it++;
    }
  }  

  // Add Cell Based Dof //
  // Get Id of Cell 
  const int id = element.getNum() * nTotVertex * nTotVertex;

  // Insert new Dof
  for(int j = 0; j < nFCell; j++){
    myDof[it] = new Dof(id, j);
    it++;
  }

  return myDof;
}

void DofManager::insertDof(Dof* d, GroupOfDof* god){
  // Try to insert Dof //
  pair<set<Dof*, DofComparator>::iterator, bool> p;
  p = dof->insert(d);
 
  // If insertion is OK (Dof 'd' didn't exist) //
  //   --> Add new Dof
  if(p.second){
    globalId->insert(pair<Dof*, int>(d, nextId));
    
    god->add(d);
    
    nextId += 1;
  }
  
  // If insertion failed (Dof 'd' already exists) //
  //   --> delete 'd' and add existing Dof
  else{
    delete d; 
    god->add(*(p.first));
  }
}
/*
void DofManager::setAsConstant(const GroupOfElement& goe, double value){
  // Get Element //
  const vector<MElement*>& element = goe.getAll();
  const unsigned int      nElement = goe.getNumber();

  // For All Element //
  for(unsigned int i = 0; i < nElement; i++){
    // Get Dof (with same key, *not* same instance)
    int   nDof;
    Dof** myDof = dofFromElement(*element[i], &nDof);

    // Look for 'Real' Dof and set Value
    for(int j = 0; j < nDof; j++){
      // Lookup in Dof set
      pair<set<Dof*, DofComparator>::iterator, bool> p 
	= dof->insert(myDof[j]);

      // If Dof doesn't exist --> Exception
      if(p.second)
	throw Exception
	  ("Dof %s don't exist: can't set value", 
	   myDof[j]->toString().c_str());

      // Else, set Dof to value
      (*(p.first))->unknown = false; // Dof is no longer an unknown
      (*(p.first))->value   = value; // And has the value 'value' 
    
      // Delete Dof
      delete myDof[j];
    }
    
    // Delete myDof
    delete[] myDof;
  }
}
*/

pair<bool, double> DofManager::getValue(const Dof& dof) const{
  map<const Dof*, double, DofComparator>::iterator end = 
    fixedDof->end();

  map<const Dof*, double, DofComparator>::iterator it = 
    fixedDof->find(&dof);

  if(it == end)
    return pair<bool, double>(false, 42);

  else
    return pair<bool, double>(true, it->second);
}

string DofManager::toString(void) const{
  stringstream s;
  map<const Dof*, int, DofComparator>::iterator i = 
    globalId->begin();
  
  map<const Dof*, int, DofComparator>::iterator end = 
    globalId->end();
  
  for(; i != end; i++)
    s << "("  << (*i).first->toString() 

      << ": " << (*i).second 
      << ")"  << endl;

  return s.str();
}
