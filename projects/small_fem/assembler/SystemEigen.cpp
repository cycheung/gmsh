#include <slepceps.h>
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

  // Alloc A //
  const size_t size = dofM->getUnfixedDofNumber();

  A = new Mat;
  MatCreate(MPI_COMM_WORLD, A);
  MatSetSizes(*A, size, size, size, size);
  MatSetType(*A, "seqaij");

  // Alloc B if needed //
  if(general){
    B = new Mat;
    MatCreate(MPI_COMM_WORLD, B);
    MatSetSizes(*B, size, size, size, size);
    MatSetType(*B, "seqaij");
  }

  // Get GroupOfDofs //
  const size_t E = fs->getSupport().getNumber();
  const vector<GroupOfDof*>& group = fs->getAllGroups();

  // Get Sparcity Pattern & PreAllocate //
  UniqueSparsity uniqueSparsity(size);
  PetscInt* nonZero = new PetscInt[size];

  // Matrix A
  for(size_t i = 0; i < size; i++)
    nonZero[i] = 0;

  for(size_t i = 0; i < E; i++)
    SystemAbstract::sparsity(nonZero, uniqueSparsity, *group[i]);

  MatSeqAIJSetPreallocation(*A, 42, nonZero);

  // Matrix B if needed (same sparsity)
  if(general)
    MatSeqAIJSetPreallocation(*B, 42, nonZero);

  delete[] nonZero;

  // Assemble SystemEigen //
  formulationPtr termA = &Formulation::weak;
  formulationPtr termB = &Formulation::weakB;

  for(size_t i = 0; i < E; i++)
    SystemAbstract::assemble(*A, i, *group[i], termA);

  if(general)
    for(size_t i = 0; i < E; i++)
      SystemAbstract::assemble(*B, i, *group[i], termB);

  MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
  if(general)
    MatAssemblyBegin(*B, MAT_FINAL_ASSEMBLY);

  MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
  if(general)
    MatAssemblyEnd(*B, MAT_FINAL_ASSEMBLY);

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
  EPSCreate(MPI_COMM_WORLD, &solver);

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
  EPSSetTolerances(solver, 1E-18, 1E6);
  EPSSetType(solver, EPSKRYLOVSCHUR);
  EPSSetWhichEigenpairs(solver, EPS_SMALLEST_MAGNITUDE);

  // Override with PETSc Database //
  EPSSetFromOptions(solver);

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
