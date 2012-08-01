#include <set>
#include <sstream>
#include "MVertex.h"
#include "GroupOfElement.h"

using namespace std;

GroupOfElement::GroupOfElement(GEntity& entity){
  // Save Entity //
  this->entity = &entity;
  
  // Get Number of Mesh Elements //
  nElement = entity.getNumMeshElements();

  // Get Elements //
  element = new vector<MElement*>(nElement);

  for(unsigned int i = 0; i < nElement; i++)
    (*element)[i] = entity.getMeshElement(i);

  // Init Other Struct //
  gov = NULL;
}

GroupOfElement::~GroupOfElement(void){
  delete element;
  
  if(gov)
    delete gov;
  /*
    GroupOfElement is *NOT* reponsible for
    deleting 'entity', niether the MElements of
    'element' !!
  */
}

int GroupOfElement::getNVertex(void) const{
  if(!gov)
    gov = new GroupOfVertex(*this);
  
  return gov->getNumber();
}

string GroupOfElement::toString(void) const{
  stringstream stream;
  
  stream << "*********************************************"    
	 << endl
	 << "* Group Of Element                          *"    
	 << endl
	 << "*********************************************" 
	 << endl << "*" 
	 << endl
	 << "* This group contains the following elements: " << endl;

  for(unsigned int i = 0; i < nElement; i++){
    stream << "*    -- ID: " 
	   << (*element)[i]->getNum() << endl;
  }
  
  stream << "*********************************************" 
	 << endl;
  
  return stream.str();
}
