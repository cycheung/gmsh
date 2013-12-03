#include "System.h"
#include "SolverMUMPS.h"

using namespace std;

System::System(const FormulationTyped<double>& formulation){
  // Check Formulation type //
  string formulationType = formulation.getType();
  string systemType      = getType();

  if(formulationType.compare(systemType) != 0)
    throw Exception("This System is of real type %s, but the Formulation is %s",
                    systemType.c_str(),
                    formulationType.c_str());

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
  formulationPtr term = &FormulationTyped<double>::weak;

  // Alloc //
  const size_t size = dofM->getUnfixedDofNumber();

  A = new SolverMatrix<double>(size, size);
  b = new SolverVector<double>(size);

  // Assemble //
  #pragma omp parallel for
  for(size_t i = 0; i < E; i++)
    SystemTyped::assemble(*A, *b, i, *group[i], term);

  // The system is assembled //
  assembled = true;
}

void System::solve(void){
  // Is the System assembled ? //
  if(!assembled)
    assemble();

  // Use SolverMUMPS //
  SolverMUMPS<double> solver;
  x = new fullVector<double>;

  solver.solve(*A, *b, *x);

  // System solved ! //
  solved = true;
}

size_t System::getNComputedSolution(void) const{
  return 1;
}

void System::getSolution(fullVector<double>& sol, size_t nSol) const{
  sol.setAsProxy(*x, 0, x->size());
}

void System::getSolution(fullVector<double>& sol) const{
  getSolution(sol, 0);
}


void System::addSolution(FEMSolution& feSol) const{
  if(!solved)
    throw Exception("System: addSolution -- System not solved");

  feSol.addCoefficients(0, 0, *fs, *dofM, *x);
}

void System::writeMatrix(std::string fileName,
                         std::string matrixName) const{
  A->writeToMatlabFile(fileName, matrixName);
}
