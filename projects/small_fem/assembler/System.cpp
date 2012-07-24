#include "System.h"
#include "Formulation.h"
#include "Exception.h"
#include "Solver.h"

using namespace std;

System::System(const Formulation& formulation){
  // Get Formulation //
  this->formulation = &formulation;

  // Get Dof Manager //
  dofM = new DofManager(formulation.fs());

  // Get DofManager Data //
  size = dofM->dofNumber();
  
  // Create System //
  A = new fullMatrix<double>(size, size);
  b = new fullVector<double>(size);
  x = new fullVector<double>(size);

  // The system is not assembled //
  isAssembled = false;
}

System::~System(void){
  delete A;
  delete b;
  delete x;
  delete dofM;
  // System is not responsible for deleting 'Formulations'
}
/*
void System::fixBC(const int physicalId, const double value){
  
  const multimap<int, Dof*>& physicals = dofM->getAllPhysicals();

  multimap<int, Dof*>::const_iterator j;

  pair<multimap<int, Dof*>::const_iterator, multimap<int, Dof*>::const_iterator>
    range;
  range = physicals.equal_range(physicalId);
  
  if((range.first == range.second) &&
     (range.first == physicals.end()))
    throw Exception("Unknown Physical");
  
  for(j = range.first; j != range.second; j++){
    // Get Dof Id
    int dofId = dofM->getGlobalId(*((*j).second));

    // We set the 'dofId'th row to zero
    for(int i = 0; i < size; i++)
      (*A)(dofId, i) = 0.0;

    // We set the 'dofId'th diagonal to one
    (*A)(dofId, dofId) = 1.0;
    
    // We also set the 'dofId'th RHS to 'value' 
    (*b)(dofId) = value;
  }
}
*/
void System::assemble(void){
  // Get GroupOfDofs //
  const std::vector<GroupOfDof*>& group = dofM->getAllGroups();
  const int E = dofM->groupNumber();

  // Assemble System //
  for(int i = 0; i < E; i++)
    assemble(*(group[i]));  
  
  // The system is assembled //
  isAssembled = true;
}

void System::solve(void){
  // Is the System assembled ? //
  if(!isAssembled)
    assemble();

  // Get dof value //
  Solver::solve(*A, *x, *b);
  /*
  // Set all Entities value //
  const vector<Dof*>* dof = &dofM->getAllDofs();
  const int N = dof->size();
  
  for(int i = 0; i < N; i++)
    dofM->getEntity(*((*dof)[i])).setValue((*x)(i));
  */
}

void System::assemble(GroupOfDof& group){
  const vector<Dof*>& dof = group.getAll();
  const int N = group.getNumber();

  for(int i = 0; i < N; i++){
    int dofI = dofM->getGlobalId(*(dof[i]));

    for(int j = 0; j < N; j++){
      int dofJ = dofM->getGlobalId(*(dof[j]));
      (*A)(dofI, dofJ) += 
	formulation->weak(i, j, group);
    }

    (*b)(dofI) += formulation->rhs(i, group);
  }
}
