#include "GroupOfDof.h"

GroupOfDof::GroupOfDof(int numberOfDof, const MElement& geoElement){
  element = &geoElement;
  
  nDof = numberOfDof;
  dof  = new std::vector<Dof*>(nDof);

  nextDof = 0;
}

GroupOfDof::~GroupOfDof(void){
  // GroupOfDofs are not responsible for
  // deleting dofs, orientations and Jacobian
  delete dof;
}

void GroupOfDof::add(Dof* dof){
  this->dof->at(nextDof) = dof;
  nextDof++;
}
