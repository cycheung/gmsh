#include <petscksp.h>
#include <petscpc.h>

#include "System.h"

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
  A      = NULL;
  b      = NULL;
  xPetsc = NULL;
  x      = NULL;

  // The system is not assembled and not solved //
  assembled = false;
  solved    = false;
}

System::~System(void){
  delete dofM;

  if(A){
    MatDestroy(A);
    delete A;
  }

  if(b){
    VecDestroy(b);
    delete b;
  }

  if(xPetsc){
    VecDestroy(xPetsc);
    delete xPetsc;
  }

  if(x){
    delete x;
  }
}

void System::assemble(void){
  // Enumerate //
  dofM->generateGlobalIdSpace();

  // Alloc //
  const size_t size = dofM->getUnfixedDofNumber();

  A      = new Mat;
  b      = new Vec;
  xPetsc = new Vec;

  // Create Matrix and Vectors //
  MatCreate(MPI_COMM_WORLD, A);
  MatSetSizes(*A, size, size, size, size);
  MatSetType(*A, "seqaij");

  VecCreate(MPI_COMM_WORLD, b);
  VecSetSizes(*b, size, size);
  VecSetType(*b, "seq");

  VecCreate(MPI_COMM_WORLD, xPetsc);
  VecSetSizes(*xPetsc, size, size);
  VecSetType(*xPetsc, "seq");

  // Get GroupOfDofs //
  const size_t E = fs->getSupport().getNumber();
  const vector<GroupOfDof*>& group = fs->getAllGroups();

  // Get Sparsity Pattern & PreAllocate //
  UniqueSparsity uniqueSparsity(size);
  PetscInt* nonZero = new PetscInt[size];

  for(size_t i = 0; i < size; i++)
    nonZero[i] = 0;

  for(size_t i = 0; i < E; i++)
    SystemAbstract::sparsity(nonZero, uniqueSparsity, *group[i]);

  MatSeqAIJSetPreallocation(*A, 42, nonZero);

  delete[] nonZero;

  // Assemble System //
  formulationPtr term = &Formulation::weak;

  for(size_t i = 0; i < E; i++)
    SystemAbstract::assemble(*A, *b, i, *group[i], term);

  MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(*b);
  MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
  VecAssemblyEnd(*b);

  /*
  PetscViewer fd;
  PetscViewerASCIIOpen(MPI_COMM_WORLD, "mat.m", &fd);
  PetscViewerSetFormat(fd, PETSC_VIEWER_ASCII_MATLAB);
  PetscObjectSetName((PetscObject)(*A), "A");
  MatView(*A, fd);
  */

  // The system is assembled //
  assembled = true;
}

void System::solve(void){
  // Is the System assembled ? //
  if(!assembled)
    assemble();

  // Build Solver //
  KSP solver;
  PC  precond;

  KSPCreate(MPI_COMM_WORLD, &solver);
  KSPSetOperators(solver, *A, *A, DIFFERENT_NONZERO_PATTERN);

  // Use MUMPS //
  KSPGetPC(solver, &precond);
  PCSetType(precond, PCLU);
  PCFactorSetMatSolverPackage(precond, MATSOLVERMUMPS);

  // Override with PETSc Database //
  KSPSetFromOptions(solver);
  PCSetFromOptions(precond);

  // Solve and Delete Solver //
  KSPSolve(solver, *b, *xPetsc);
  KSPDestroy(&solver);

  // Get Solution //
  double* solution;
  VecGetArray(*xPetsc, &solution);

  x = new fullVector<double>(solution, dofM->getUnfixedDofNumber());

  // System solved ! //
  solved = true;
}
