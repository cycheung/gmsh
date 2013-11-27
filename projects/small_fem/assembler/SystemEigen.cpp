#include "slepceps.h"
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

  // Init //
  A           = NULL;
  B           = NULL;
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

  if(A){
    MatDestroy(A);
    delete A;
  }

  if(B){
    MatDestroy(B);
    delete B;
  }

  delete dofM;
}

void SystemEigen::
setNumberOfEigenValues(size_t nEigenValues){
  const size_t nDof = dofM->getUnfixedDofNumber();

  if(nEigenValues > nDof)
    throw
      Exception
      ("I can't compute more Eigenvalues (%d) than the number of unknowns (%d)",
       nEigenValues, nDof);

  else
    this->nEigenValues = nEigenValues;
}

void SystemEigen::assemble(void){
  // Enumerate //
  dofM->generateGlobalIdSpace();

  // Get GroupOfDofs //
  const size_t E = fs->getSupport().getNumber();
  const vector<GroupOfDof*>& group = fs->getAllGroups();

  // Get Formulation Terms //
  formulationPtr termA = &Formulation::weak;
  formulationPtr termB = &Formulation::weakB;

  // Alloc Temp Sparse Matrices (not with PETSc) //
  const size_t size = dofM->getUnfixedDofNumber();

  SolverVector tmpRHS(size);
  SolverMatrix tmpA(size, size);
  SolverMatrix tmpB(size, size);

  // Assemble Systems (tmpA and tmpB) //
  #pragma omp parallel for
  for(size_t i = 0; i < E; i++)
    SystemAbstract::assemble(tmpA, tmpRHS, i, *group[i], termA);

  if(general)
    #pragma omp parallel for
    for(size_t i = 0; i < E; i++)
      SystemAbstract::assemble(tmpB, tmpRHS, i, *group[i], termB);

  // Copy tmpA into Assembled PETSc matrix //
  // Data
  vector<int>    row;
  vector<int>    col;
  vector<double> value;
  int            nNZ;

  // Serialize (CStyle) tmpA & Copy
  nNZ = tmpA.serializeCStyle(row, col, value);
  A   = new Mat;

  MatCreateSeqAIJFromTriple(MPI_COMM_SELF, size, size,
                            row.data(), col.data(), value.data(),
                            A, nNZ, PETSC_FALSE);

  // Copy tmpB (CStyle) into Assembled PETSc matrix (if needed) //
  if(general){
    nNZ = tmpB.serializeCStyle(row, col, value);
    B   = new Mat;

    MatCreateSeqAIJFromTriple(MPI_COMM_SELF, size, size,
                              row.data(), col.data(), value.data(),
                              B, nNZ, PETSC_FALSE);
  }

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

  // Build Solver //
  EPS solver;
  EPSCreate(MPI_COMM_SELF, &solver);

  if(general)
    EPSSetOperators(solver, *A, *B);
  else
    EPSSetOperators(solver, *A, NULL);

  if(general)
    EPSSetProblemType(solver, EPS_GNHEP);
  else
    EPSSetProblemType(solver, EPS_NHEP);

  // Set Options //
  EPSSetDimensions(solver, nEigenValues, PETSC_DECIDE, PETSC_DECIDE);
  EPSSetTolerances(solver, 1E-12, 1E6);
  EPSSetWhichEigenpairs(solver, EPS_SMALLEST_MAGNITUDE);

  // Use Krylov Schur //
  EPSSetType(solver, EPSKRYLOVSCHUR);
  /*
  // Use Generalized Davidson Solver and LU (MUMPS) preconditioning //
  KSP linSolver;
  PC  precond;
  ST  specT;

  EPSSetType(solver, "gd");

  EPSGetST(solver, &specT);
  STSetType(specT, "precond");
  STGetKSP(specT, &linSolver);

  KSPSetType(linSolver, "preonly");
  KSPGetPC(linSolver, &precond);
  PCSetType(precond, "lu");
  PCFactorSetMatSolverPackage(precond, "mumps");
  */

  // Override with PETSc Database //
  EPSSetFromOptions(solver);
  //STSetFromOptions(specT);

  // Solve //
  EPSSolve(solver);

  // Get Solution //
  const size_t size = dofM->getUnfixedDofNumber();

  PetscScalar  lambdaReal;
  PetscScalar  lambdaImag;
  PetscScalar* xReal;
  PetscScalar* xImag;
  Vec          xRealPetsc;
  Vec          xImagPetsc;

  MatGetVecs(*A, PETSC_NULL, &xRealPetsc);
  MatGetVecs(*A, PETSC_NULL, &xImagPetsc);

  EPSGetConverged(solver, &nEigenValues);

  eigenValue  = new vector<complex<double> >(nEigenValues);
  eigenVector = new vector<fullVector<complex<double> > >(nEigenValues);

  for(PetscInt i = 0; i < nEigenValues; i++){
    EPSGetEigenpair(solver, i,
                    &lambdaReal, &lambdaImag,
                    xRealPetsc, xImagPetsc);

    VecGetArray(xRealPetsc, &xReal);
    VecGetArray(xImagPetsc, &xImag);

    (*eigenVector)[i].resize(size);
    for(size_t j = 0; j < size; j++)
      (*eigenVector)[i](j) = complex<double>(xReal[j], xImag[j]);

    (*eigenValue)[i] = complex<double>(lambdaReal, lambdaImag);
  }

  VecDestroy(&xRealPetsc);
  VecDestroy(&xImagPetsc);
  EPSDestroy(&solver);

  // System solved ! //
  solved = true;
}

void SystemEigen::addSolution(FEMSolution& feSol) const{
  if(!solved)
    throw Exception("System: addSolution -- System not solved");

  for(int i = 0; i < nEigenValues; i++)
    feSol.addCoefficients
      (i * 2, 0, *fs, *dofM, (*eigenVector)[i]);
}
