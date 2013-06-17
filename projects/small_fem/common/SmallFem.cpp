#include <mpi.h>
#include <petsc.h>
#include <slepc.h>

#include "BasisFactory.h"
#include "Context.h"

#include "SmallFem.h"

bool SmallFem::initOne = false;
bool SmallFem::finaOne = false;

SmallFem::SmallFem(void){
}

SmallFem::~SmallFem(void){
}

void SmallFem::Initialize(int argc, char** argv){
  // Initialize only once

  if(!initOne){
    // Call MPI, PETSc and SLEPc
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

    // Gmsh Instance
    CTX::instance();

    initOne = true;
  }
}

void SmallFem::Finalize(void){
  // Finalize only once and if Initialize first

  if(!finaOne && initOne){
    // Clear Gmsh Instance
    CTX* tmp = CTX::instance();

    if(tmp)
      delete tmp;

    // Clear Gmsh BasisFactory
    BasisFactory::clearAll();

    // Finalize MPI, PETSc and SLEPc
    PetscFinalize();
    SlepcFinalize();
    MPI_Finalize();

    finaOne = true;
  }
}
