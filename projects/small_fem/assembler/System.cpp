#include "System.h"
#include "Formulation.h"
#include "Exception.h"
#include "Solver.h"

using namespace std;

System::System(const Formulation& formulation){
  // Get Formulation //
  this->formulation = &formulation;
  this->fs          = &(formulation.fs());

  // Get Dof Manager //
  dofM = new DofManager(*fs);

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

  // The system is assembled //
  isAssembled = true;  
}

void System::fixBC(const GroupOfElement& goe, double value){
  const vector<MElement*>&  element = goe.getAll();
  unsigned int             nElement = goe.getNumber();
  
  for(unsigned int i = 0; i < nElement; i++){
    const vector<Dof*>   dof = fs->getKeys(*element[i]);
    const unsigned int nDof  = dof.size();
    
    for(unsigned int j = 0; j < nDof; j++){
      dofM->fixValue(*dof[j], value);

      delete dof[j];
    }
  }
}

void System::solve(void){
  // Is the System assembled ? //
  if(!isAssembled)
    assemble();

  // Get dof value //
  Solver::solve(*A, *x, *b);
}
			   
void System::assemble(GroupOfDof& group){
  const vector<Dof*>& dof = group.getAll();
  const int N = group.getNumber();

  for(int i = 0; i < N; i++){
    pair<bool, double> fixed = dofM->getValue(*(dof[i]));
    int dofI = dofM->getGlobalId(*(dof[i]));
    
    if(fixed.first){
      // If fixed Dof
      (*A)(dofI, dofI) = 1;
      (*b)(dofI)       = fixed.second;
    }
       
    else{
      // If unknown Dof
      for(int j = 0; j < N; j++){
	int dofJ = dofM->getGlobalId(*(dof[j]));
	
	(*A)(dofI, dofJ) += 
	  formulation->weak(i, j, group);
      }
      
      (*b)(dofI) += formulation->rhs(i, group);
    }
  } 
}

