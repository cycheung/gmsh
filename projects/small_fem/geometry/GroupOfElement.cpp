#include <set>
#include <sstream>
#include "MVertex.h"
#include "GroupOfElement.h"

using namespace std;

unsigned int GroupOfElement::nextId = 0;

GroupOfElement::GroupOfElement(GEntity& entity, Mesh& mesh){
  // Set ID //
  id = nextId;
  nextId++;
  
  // Save Entity & Mesh//
  this->entity = &entity;
  this->mesh   = &mesh;

  // Get Number of Mesh Elements //
  nElement = entity.getNumMeshElements();

  // Get Elements //
  element = new vector<MElement*>(nElement);

  for(unsigned int i = 0; i < nElement; i++)
    (*element)[i] = entity.getMeshElement(i);

  // Init Other Struct //
  gov = NULL;
  goe = NULL;
}

GroupOfElement::~GroupOfElement(void){
  delete element;

  if(goe)
    delete goe;
}

GroupOfVertex& GroupOfElement::getGroupOfVertex(void){
  if(!gov)
    gov = &(mesh->getGroupOfVertex(*this));
  
  return *gov;
}

GroupOfEdge& GroupOfElement::getGroupOfEdge(void){
  if(!goe)
    goe = &(mesh->getGroupOfEdge(*this));
  
  return *goe;
}

string GroupOfElement::toString(void) const{
  stringstream stream;
  
  stream << "*********************************************"    
	 << endl
	 << "* Group Of Element #" << id    
	 << endl
	 << "*********************************************" 
	 << endl << "*" 
	 << endl
	 << "* This group contains the following elements: " << endl;

  for(unsigned int i = 0; i < nElement; i++)
    stream << "*    -- ID: " 
	   << mesh->getGlobalId(*(*element)[i]) << endl;
 
  stream << "*********************************************" 
	 << endl;
  
  return stream.str();
}
