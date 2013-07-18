#include "SystemEigen.h"

using namespace std;

SystemEigen::SystemEigen(const Formulation& formulation){
  // Get Formulation //
  this->formulation = &formulation;
  this->fs          = &(formulation.fs());

  // Get Dof Manager //
  dofM = new DofManager();
  dofM->addToDofManager(fs->getAllGroups());

  // Is the Problem a General EigenValue Problem ? //
  general = formulation.isGeneral();

  // Create SystemEigen //

  // Init
  linSysA = NULL;
  linSysB = NULL;

  // eSys will be created at solving point
  eSys        = NULL;
  eigenValue  = NULL;
  eigenVector = NULL;

  // The SystemEigen is not assembled and not solved//
  nEigenValues = 0;
  assembled    = false;
  solved       = false;
}

SystemEigen::~SystemEigen(void){
  if(eigenVector)
    delete eigenVector;

  if(eigenValue)
    delete eigenValue;

  if(eSys)
    delete eSys;

  if(linSysA)
    delete linSysA;

  if(general)
    delete linSysB;

  delete dofM;
  // SystemEigen is not responsible for deleting 'Formulations'
}

void SystemEigen::
setNumberOfEigenValues(size_t nEigenValues){
  const size_t nDof = dofM->getDofNumber();

  if(nEigenValues > nDof)
    throw
      Exception("I cannot compute more Eigenvalues (%d) than the number of unknowns (%d)",
                nEigenValues, nDof);

  else
    this->nEigenValues = nEigenValues;
}

void SystemEigen::assemble(void){
  const size_t size = dofM->getDofNumber();

  // Enumerate //
  dofM->generateGlobalIdSpace();

  // Init System //
  // Linear System A
  linSysA = new linearSystemPETSc<double>();
  linSysA->allocate(size);

  // Linear System B
  if(general){
    linSysB = new linearSystemPETSc<double>();
    linSysB->allocate(size);
  }

  // Get GroupOfDofs //
  const size_t E = fs->getSupport().getNumber();
  const vector<GroupOfDof*>& group = fs->getAllGroups();

  // Get Sparcity Pattern & PreAllocate//
  // linSysA
  for(size_t i = 0; i < E; i++)
    SystemAbstract::sparsity(*linSysA, *group[i]);

  linSysA->preAllocateEntries();

  // linSysB
  if(general){
    for(size_t i = 0; i < E; i++)
      SystemAbstract::sparsity(*linSysB, *group[i]);

    linSysB->preAllocateEntries();
  }

  // Assemble SystemEigen //
  formulationPtr termA = &Formulation::weak;
  formulationPtr termB = &Formulation::weakB;

  for(size_t i = 0; i < E; i++)
    SystemAbstract::assemble(*linSysA, i, *group[i], termA);

  if(general)
    for(size_t i = 0; i < E; i++)
      SystemAbstract::assemble(*linSysB, i, *group[i], termB);

  // The SystemEigen is assembled //
  assembled = true;
}

void SystemEigen::solve(void){
  // Check nEigenValues
  if(!nEigenValues)
    throw
      Exception("The number of eigenvalues to compute is zero");

  // Is the SystemEigen assembled ? //
  if(!assembled)
    assemble();

  // Solve //
  eSys = new EigenSolver(linSysA, linSysB, false);
  eSys->solve(nEigenValues, "smallest");

  // Get Solution //
  nEigenValues = eSys->getNumEigenValues();

  eigenValue  = new vector<complex<double> >(nEigenValues);
  eigenVector = new vector<vector<complex<double> > >(nEigenValues);

  for(size_t i = 0; i < nEigenValues; i++)
    (*eigenValue)[i] = eSys->getEigenValue(i);

  for(size_t i = 0; i < nEigenValues; i++)
    (*eigenVector)[i] = eSys->getEigenVector(i);

  // System solved ! //
  solved = true;
}
