#include "System.h"
#include "SolverMUMPS.h"

using namespace std;

System::System(const Formulation<double>& formulation){
  // Get Formulation //
  this->formulation = &formulation;
  this->fs          = &(formulation.fs());

  // Get Dof Manager //
  dofM = new DofManager<double>();
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
  formulationPtr term = &Formulation<double>::weak;

  // Alloc //
  const size_t size = dofM->getUnfixedDofNumber();

  A = new SolverMatrix<double>(size, size);
  b = new SolverVector<double>(size);

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

void System::getSolution(std::map<Dof, double>& sol, size_t nSol) const{
  // Get All Dofs
  const vector<Dof> dof = fs->getAllDofs();
  const size_t     nDof = dof.size();

  // Fill Map
  for(size_t i = 0; i < nDof; i++)
    sol.insert(pair<Dof, double>(dof[i], (*x)(dofM->getGlobalId(dof[i]))));
}

void System::getSolution(std::map<Dof, double>& sol) const{
  getSolution(sol, 0);
}

void System::getSolution(FEMSolution<double>& feSol) const{
  if(!solved)
    throw Exception("System: addSolution -- System not solved");

  feSol.addCoefficients(0, 0, *fs, *dofM, *x);
}

void System::writeMatrix(std::string fileName,
                         std::string matrixName) const{
  A->writeToMatlabFile(fileName, matrixName);
}
