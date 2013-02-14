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
  const unsigned int E = fs->getSupport().getNumber();
  const vector<pair<const MElement*, ElementData> >&
    element = fs->getSupport().getAll();

  // Get Sparsity Pattern & PreAllocate//
  for(unsigned int i = 0; i < E; i++)
    SystemAbstract::sparsity(*linSys,
                             element[i].second.getGroupOfDof());

  linSys->preAllocateEntries();

  // Assemble System //
  formulationPtr term = &Formulation::weak;

  for(unsigned int i = 0; i < E; i++)
    SystemAbstract::assemble(*linSys,
                             element[i].second.getGroupOfDof(),
                             term);

  // The system is assembled //
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
