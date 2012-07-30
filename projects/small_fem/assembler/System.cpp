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

void System::assemble(void){
  // Get GroupOfDofs //
  const std::vector<GroupOfDof*>& group = dofM->getAllGroups();
  const int E = dofM->groupNumber();

  // Assemble System //
  for(int i = 0; i < E; i++)
    assemble(*(group[i]));  
  
  // Set BC //
  for(int i = 0; i < size; i++){
    (*A)(0, i) = 0.0;
    (*A)(1, i) = 0.0;
    (*A)(4, i) = 0.0;
    (*A)(5, i) = 0.0;
    (*A)(6, i) = 0.0;
    (*A)(8, i) = 0.0;
  }

  (*A)(0, 0) = 1.0;
  (*A)(1, 1) = 1.0;
  (*A)(4, 4) = 1.0;
  (*A)(5, 5) = 1.0;
  (*A)(6, 6) = 1.0;
  (*A)(8, 8) = 1.0;  

  (*b)(0) = -1.0;
  (*b)(1) = -1.0;
  (*b)(4) =  2.0;
  (*b)(5) =  2.0;
  (*b)(6) = -1.0;
  (*b)(8) =  2.0;    

  // The system is assembled //
  isAssembled = true;  
}

void System::fixBC(const GroupOfElement& goe, double value){
  //dofM->setAsConstant(goe, value);
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

