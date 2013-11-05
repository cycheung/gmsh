#include "System.h"
#include "SolverMUMPS.h"

using namespace std;

System::System(void){
  // Do nothing //
  // Just for inheritance //
}

System::System(const Formulation& formulation){
  // Get Formulation //
  this->formulation = &formulation;
  this->fs          = &(formulation.fs());

  // Get Dof Manager //
  dofM = new DofManager();
  dofM->addToDofManager(fs->getAllGroups());

  // Init //
  A = NULL;
  b = NULL;
  x = NULL;

  // The system is not assembled and not solved //
  assembled = false;
  solved    = false;
}

System::~System(void){
  delete dofM;

  if(A)
    delete A;

  if(b)
    delete b;

  if(x)
    delete x;
}

void System::assemble(void){
  // Enumerate //
  dofM->generateGlobalIdSpace();

  // Get GroupOfDofs //
  const size_t E = fs->getSupport().getNumber();
  const vector<GroupOfDof*>& group = fs->getAllGroups();

  // Get Formulation Term //
  formulationPtr term = &Formulation::weak;

  // Alloc //
  const size_t size = dofM->getUnfixedDofNumber();

  A = new SolverMatrix(size, size);
  b = new SolverVector(size);

  // Assemble //
  #pragma omp parallel for
  for(size_t i = 0; i < E; i++)
    SystemAbstract::assemble(*A, *b, i, *group[i], term);

  // The system is assembled //
  assembled = true;
}

void System::solve(void){
  // Is the System assembled ? //
  if(!assembled)
    assemble();

  // Use SolverMUMPS //
  SolverMUMPS solver;
  x = new fullVector<double>;

  solver.solve(*A, *b, *x);

  // System solved ! //
  solved = true;
}
