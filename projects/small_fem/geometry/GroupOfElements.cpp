#include <sstream>
#include "GroupOfElements.h"

using namespace std;

GroupOfElements::GroupOfElements(GEntity& entity, int id){
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

GroupOfElements::~GroupOfElements(void){
  delete element;
  
  /*
    GroupOfElements is *NOT* reponsible for
    deleting 'entity', niether the MElements of
    'element' !!
  */
}

string GroupOfElements::toString(void) const{
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
