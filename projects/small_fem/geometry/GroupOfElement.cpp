#include <sstream>
#include "GroupOfElement.h"

using namespace std;

GroupOfElement::GroupOfElement(GEntity& entity, int id){
  // Save Entity //
  this->id     = id;
  this->entity = &entity;
  
  // Get Number of Mesh Elements //
  nElement = entity.getNumMeshElements();

  // Get Elements //
  element = new vector<MElement*>(nElement);

  for(unsigned int i = 0; i < nElement; i++)
    (*element)[i] = entity.getMeshElement(i);
}

GroupOfElement::~GroupOfElement(void){
  delete element;
  
  /*
    GroupOfElement is *NOT* reponsible for
    deleting 'entity', niether the MElements of
    'element' !!
  */
}

string GroupOfElement::toString(void) const{
  stringstream stream;
  
  stream << "*********************************************"    
	 << endl
	 << "* Groups Of Elements number " << id    
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
