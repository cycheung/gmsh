#include "System.h"

using namespace std;

System::System(const Formulation& formulation){
  // Get Formulation //
  this->formulation = &formulation;
  this->fs          = &(formulation.fs());

  // Get Dof Manager //
  dofM = new DofManager();
  dofM->addToDofManager(fs->getAllGroups());

  // Init //
  x      = NULL;
  linSys = NULL;

  // The system is not assembled and not solved //
  assembled = false;
  solved    = false;
}

System::~System(void){
  delete x;
  delete dofM;
  if(linSys)
    delete linSys;
  // System is not responsible for deleting 'Formulations'
}

void System::assemble(void){
  // Enumerate //
  dofM->generateGlobalIdSpace();

  // Init System //
  x      = new fullVector<double>(dofM->getDofNumber());
  linSys = new linearSystemPETSc<double>();
  linSys->allocate(dofM->getDofNumber());

  // Get GroupOfDofs //
  const unsigned int E = fs->getSupport().getNumber();
  const vector<GroupOfDof*>& group = fs->getAllGroups();

  // Get Sparsity Pattern & PreAllocate//
  for(unsigned int i = 0; i < E; i++)
    SystemAbstract::sparsity(*linSys, *group[i]);

  linSys->preAllocateEntries();

  // Assemble System //
  formulationPtr term = &Formulation::weak;

  for(unsigned int i = 0; i < E; i++)
    SystemAbstract::assemble(*linSys, i, *group[i], term);

  // The system is assembled //
  assembled = true;
}

void System::solve(void){
  // Is the System assembled ? //
  if(!assembled)
    assemble();

  // Solve //
  linSys->systemSolve();
  //linSys->print();
  // Write Sol
  const unsigned int size = dofM->getDofNumber();
  double xi;

  for(unsigned int i = 0; i < size; i++){
    linSys->getFromSolution(i, xi);
    (*x)(i) = xi;
  }

  // System solved ! //
  solved = true;
}
