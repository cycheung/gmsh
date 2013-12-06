///////////////////////////////////////////
// Templates Implementations for System: //
// Inclusion compilation model           //
//                                       //
// Damn you gcc: we want 'export' !      //
///////////////////////////////////////////

#include "SolverMUMPS.h"

template<typename scalar>
System<scalar>::System(const Formulation<scalar>& formulation){
  // Get Formulation //
  this->formulation = &formulation;
  this->fs          = &(formulation.fs());

  // Get Dof Manager //
  this->dofM = new DofManager<scalar>();
  this->dofM->addToDofManager(this->fs->getAllGroups());

  // Init //
  A = NULL;
  b = NULL;
  x = NULL;

  // The system is not assembled and not solved //
  this->assembled = false;
  this->solved    = false;
}

template<typename scalar>
System<scalar>::~System(void){
  delete this->dofM;

  if(A)
    delete A;

  if(b)
    delete b;

  if(x)
    delete x;
}

template<typename scalar>
void System<scalar>::addBorderTerm(const Formulation<scalar>& formulation){
  // Get the FunctionSpace of the given formulation
  const FunctionSpace& fs = formulation.fs();

  // Get GroupOfDofs //
  const size_t                        E = fs.getSupport().getNumber();
  const std::vector<GroupOfDof*>& group = fs.getAllGroups();

  // Get Formulation Term //
  typename SystemAbstract<scalar>::formulationPtr term =
    &Formulation<scalar>::weak;

  // Assemble //
  #pragma omp parallel for
  for(size_t i = 0; i < E; i++)
    SystemAbstract<scalar>::
      assemble(*A, *b, i, *group[i], term, formulation);
}

template<typename scalar>
void System<scalar>::assemble(void){
  // Enumerate //
  this->dofM->generateGlobalIdSpace();

  // Get GroupOfDofs //
  const size_t                        E = this->fs->getSupport().getNumber();
  const std::vector<GroupOfDof*>& group = this->fs->getAllGroups();

  // Get Formulation Term //
  typename SystemAbstract<scalar>::formulationPtr term =
    &Formulation<scalar>::weak;

  // Alloc //
  const size_t size = this->dofM->getUnfixedDofNumber();

  A = new SolverMatrix<scalar>(size, size);
  b = new SolverVector<scalar>(size);

  // Assemble //
  #pragma omp parallel for
  for(size_t i = 0; i < E; i++)
    SystemAbstract<scalar>::
      assemble(*A, *b, i, *group[i], term, *this->formulation);

  // The system is assembled //
  this->assembled = true;
}

template<typename scalar>
void System<scalar>::solve(void){
  // Is the System assembled ? //
  if(!this->assembled)
    assemble();

  // Use SolverMUMPS //
  SolverMUMPS<scalar> solver;
  x = new fullVector<scalar>;

  solver.solve(*A, *b, *x);

  // System solved ! //
  this->solved = true;
}

template<typename scalar>
size_t System<scalar>::getNComputedSolution(void) const{
  return 1;
}

template<typename scalar>
void System<scalar>::getSolution(fullVector<scalar>& sol, size_t nSol) const{
  sol.setAsProxy(*x, 0, x->size());
}

template<typename scalar>
void System<scalar>::getSolution(fullVector<scalar>& sol) const{
  getSolution(sol, 0);
}

template<typename scalar>
void System<scalar>::getSolution(std::map<Dof, scalar>& sol, size_t nSol) const{
  // Get All Dofs
  const std::vector<Dof> dof = this->fs->getAllDofs();
  const size_t          nDof = dof.size();

  // Fill Map
  for(size_t i = 0; i < nDof; i++)
    sol.insert(std::pair<Dof, scalar>
               (dof[i], (*x)(this->dofM->getGlobalId(dof[i]))));
}

template<typename scalar>
void System<scalar>::getSolution(std::map<Dof, scalar>& sol) const{
  getSolution(sol, 0);
}

template<typename scalar>
void System<scalar>::getSolution(FEMSolution<scalar>& feSol) const{
  if(!this->solved)
    throw Exception("System: addSolution -- System not solved");

  feSol.addCoefficients(0, 0, *this->fs, *this->dofM, *x);
}

template<typename scalar>
void System<scalar>::writeMatrix(std::string fileName,
                                 std::string matrixName) const{
  A->writeToMatlabFile(fileName, matrixName);
}
