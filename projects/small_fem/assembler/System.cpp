#include "System.h"
#include "Formulation.h"
#include "Exception.h"
#include "Solver.h"

#include <cstdio>

using namespace std;

System::System(const std::vector<Element*>& elements,
	       const Formulation& formulation){
  // Get Formulation //
  this->formulation = &formulation;

  // Get Dof Manager //
  dofM = new DofManager(elements);

  // Get DofManager Data //
  size = dofM->dofNumber();
  
  const std::vector<GroupOfDof*>& group = dofM->getAllGroups();
  const int E = dofM->groupNumber();

  // Create System //
  A = new Matrix(size, size);
  n = new Vector<double>(size);

  A->allToZero();
  n->allToZero();

  // Assemble System //
  for(int i = 0; i < E; i++)
    assemble(*(group[i]));
}

System::~System(void){
  delete A;
  delete n;
  delete dofM;
  // System is not responsible for deleting 'Formulations'
}

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
    (*n)(dofId) = value;
  }
}

void System::solve(void){
  // Get dof value //
  Solver::solve(*A, *n);

  // Set all Entities value //
  const vector<Dof*>* dof = &dofM->getAllDofs();
  const int N = dof->size();
  
  for(int i = 0; i < N; i++)
    dofM->getEntity(*((*dof)[i])).setValue((*n)(i));
}

void System::assemble(GroupOfDof& group){
  const vector<Dof*>& dof = group.getAllDofs();
  const int N = group.dofNumber();

  for(int i = 0; i < N; i++){
    int dofI = dofM->getGlobalId(*(dof[i]));

    for(int j = 0; j < N; j++){
      int dofJ = dofM->getGlobalId(*(dof[j]));
      (*A)(dofI, dofJ) += 
	formulation->weak(i, j, group);
    }

    (*n)(dofI) += formulation->rhs(i, group);
  }
}
