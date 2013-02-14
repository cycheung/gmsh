#include "EigenSystem.h"

using namespace std;

EigenSystem::EigenSystem(const Formulation& formulation){
  // Get Formulation //
  this->formulation = &formulation;
  this->fs          = &(formulation.fs());

  // Get Dof Manager //
  dofM = new DofManager();
  dofM->addToGlobalIdSpace(fs->getAllGroups());

  // Is the Problem a General EigenValue Problem ? //
  general = formulation.isGeneral();

  // Create EigenSystem //
  // Get Size
  const unsigned int size = fs->dofNumber();

  // Linear System A
  linSysA = new linearSystemPETSc<double>();
  linSysA->allocate(size);

  // Linear System B
  if(general){
    linSysB = new linearSystemPETSc<double>();
    linSysB->allocate(size);
  }

  else{
    linSysB = NULL;
  }

  // eSys will be created at solving point
  eSys        = NULL;
  eigenValue  = NULL;
  eigenVector = NULL;

  // The EigenSystem is not assembled and not solved//
  nEigenValues = 0;
  assembled    = false;
  solved       = false;
}

EigenSystem::~EigenSystem(void){
  if(eigenVector)
    delete eigenVector;

  if(eigenValue)
    delete eigenValue;

  if(eSys)
    delete eSys;

  delete linSysA;

  if(general)
    delete linSysB;

  delete dofM;
  // EigenSystem is not responsible for deleting 'Formulations'
}

void EigenSystem::
setNumberOfEigenValues(unsigned int nEigenValues){
  const unsigned int nDof = fs->dofNumber();

  if(nEigenValues > nDof)
    throw
      Exception("I cannot compute more Eigenvalues (%d) than the number of unknowns (%d)",
		nEigenValues, nDof);

  else
    this->nEigenValues = nEigenValues;
}

void EigenSystem::assemble(void){
  // Get GroupOfDofs //
  const vector<GroupOfDof*>& group = fs->getAllGroups();
  const unsigned int E = fs->groupNumber();

  // Get Sparcity Pattern & PreAllocate//
  // linSysA
  for(unsigned int i = 0; i < E; i++)
    AbstractSystem::sparsity(*linSysA, *(group[i]));

  linSysA->preAllocateEntries();

  // linSysB
  if(general){
    for(unsigned int i = 0; i < E; i++)
      AbstractSystem::sparsity(*linSysB, *(group[i]));

    linSysB->preAllocateEntries();
  }

  // Assemble EigenSystem //
  formulationPtr termA = &Formulation::weak;
  formulationPtr termB = &Formulation::weakB;

  for(unsigned int i = 0; i < E; i++)
    AbstractSystem::assemble(*linSysA, *(group[i]), termA);

  if(general)
    for(unsigned int i = 0; i < E; i++)
      AbstractSystem::assemble(*linSysB, *(group[i]), termB);

  // The EigenSystem is assembled //
  assembled = true;
}

void EigenSystem::solve(void){
  // Check nEigenValues
  if(!nEigenValues)
    throw
      Exception("The number of eigenvalues to compute is zero");

  // Is the EigenSystem assembled ? //
  if(!assembled)
    assemble();

  // Solve //
  eSys = new EigenSolver(linSysA, linSysB, false);
  eSys->solve(nEigenValues, "smallest");

  // Get Solution //
  nEigenValues = eSys->getNumEigenValues();

  eigenValue  = new vector<complex<double> >(nEigenValues);
  eigenVector = new vector<vector<complex<double> > >(nEigenValues);

  for(unsigned int i = 0; i < nEigenValues; i++)
    (*eigenValue)[i] = eSys->getEigenValue(i);

  for(unsigned int i = 0; i < nEigenValues; i++)
    (*eigenVector)[i] = eSys->getEigenVector(i);

  // System solved ! //
  solved = true;
}
