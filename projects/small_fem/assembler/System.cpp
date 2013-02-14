#include "System.h"

using namespace std;

System::System(const Formulation& formulation){
  // Get Formulation //
  this->formulation = &formulation;
  this->fs          = &(formulation.fs());

  // Get Dof Manager //
  dofM = new DofManager();
  dofM->addToGlobalIdSpace(fs->getAllGroups());

  // Create System //
  const unsigned int size = fs->dofNumber();

  x      = new fullVector<double>(size);
  linSys = new linearSystemPETSc<double>();
  linSys->allocate(size);

  // The system is not assembled and not solved //
  assembled = false;
  solved    = false;
}

System::~System(void){
  delete x;
  delete linSys;
  delete dofM;
  // System is not responsible for deleting 'Formulations'
}

void System::assemble(void){
  // Get GroupOfDofs //
  const vector<GroupOfDof*>& group = fs->getAllGroups();
  const int E = fs->groupNumber();

  // Set to put Fixed Dof only ones
  // (cannot use both 'setValue' and 'addValue' in PETSc)
  fixedOnes = new set<const Dof*, DofComparator>();

  // Get Sparsity Pattern & PreAllocate//
  for(int i = 0; i < E; i++)
    sparsity(*(group[i]));

  linSys->preAllocateEntries();

  // Assemble System //
  for(int i = 0; i < E; i++)
    assemble(*(group[i]));

  // The system is assembled //
  delete fixedOnes;
  assembled = true;
}

void System::solve(void){
  // Is the System assembled ? //
  if(!assembled)
    assemble();

  // Solve //
  linSys->systemSolve();

  // Write Sol
  const unsigned int size = fs->dofNumber();
  double xi;

  for(unsigned int i = 0; i < size; i++){
    linSys->getFromSolution(i, xi);
    (*x)(i) = xi;
  }

  // System solved ! //
  solved = true;
}

void System::assemble(GroupOfDof& group){
  const vector<const Dof*>& dof = group.getAll();
  const int N = group.getNumber();

  for(int i = 0; i < N; i++){
    pair<bool, double> fixed = dofM->getValue(*(dof[i]));
    int dofI = dofM->getGlobalId(*(dof[i]));

    if(fixed.first){
      // If fixed Dof
      pair<
	set<const Dof*, DofComparator>::iterator,
	bool> ones = fixedOnes->insert(dof[i]);

      if(ones.second){
	linSys->addToMatrix(dofI, dofI, 1);
	linSys->addToRightHandSide(dofI, fixed.second);
      }
    }

    else{
      // If unknown Dof
      for(int j = 0; j < N; j++){
	int dofJ = dofM->getGlobalId(*(dof[j]));

	linSys->addToMatrix(dofI, dofJ,
			    formulation->weak(i, j, group));
      }

      linSys->addToRightHandSide(dofI,
				 formulation->rhs(i, group));
    }
  }
}

void System::sparsity(GroupOfDof& group){
  const vector<const Dof*>& dof = group.getAll();
  const int N = group.getNumber();

  for(int i = 0; i < N; i++){
    pair<bool, double> fixed = dofM->getValue(*(dof[i]));
    int dofI = dofM->getGlobalId(*(dof[i]));

    if(fixed.first){
      // If fixed Dof
      linSys->insertInSparsityPattern(dofI, dofI);
    }

    else{
      // If unknown Dof
      for(int j = 0; j < N; j++){
	int dofJ = dofM->getGlobalId(*(dof[j]));

	linSys->insertInSparsityPattern(dofI, dofJ);
      }
    }
  }
}
