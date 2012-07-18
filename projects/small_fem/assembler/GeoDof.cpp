#include "GeoDof.h"

GeoDof::GeoDof(int numberOfDof, const MElement& geoElement){
  element = const_cast<MElement*>(&geoElement);
  
  nDof = numberOfDof;
  dof  = new std::vector<Dof*>(nDof);

  nextDof = 0;
}

GeoDof::~GeoDof(void){
  // GeoDofs are not responsible for
  // deleting dofs, orientations and Jacobian
  delete dof;
}

void GeoDof::add(Dof* dof){
  this->dof->at(nextDof) = dof;
  nextDof++;
}
