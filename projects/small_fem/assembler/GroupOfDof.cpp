#include "GroupOfDof.h"

GroupOfDof::GroupOfDof(int numberOfDof, int groupId){
  id = groupId;
  
  nextDof = 0;
  nDof = numberOfDof;
  dof = new std::vector<Dof*>(nDof);
  jac = NULL;
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
