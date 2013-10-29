//#include <mpi.h>
//#include <petsc.h>
//#include <slepc.h>

#include <vector>
#include <string>

#include "BasisFactory.h"
#include "Context.h"

#include "SmallFem.h"

bool     SmallFem::initOne = false;
bool     SmallFem::finaOne = false;
Options* SmallFem::option  = NULL;

SmallFem::SmallFem(void){
}

SmallFem::~SmallFem(void){
  if(option)
    delete option;
}

void SmallFem::Initialize(int argc, char** argv){
  // Initialize only once

  if(!initOne){
    // Call MPI
    //MPI_Init(&argc, &argv);

    // Get Options (Remove command name)
    option = new Options(argc - 1, argv + 1);

    // PETSc And SLEPc
    /*
    std::vector<std::string> argPetsc = option->getValue("-solver");
    int argPetscSize = argPetsc.size() + 1;

    char** argCStylePetsc = new char*[argPetscSize];

    argCStylePetsc[0] = argv[0];
    Options::cStyle(argPetsc, argCStylePetsc, 1);

    PetscInitialize(&argPetscSize, &argCStylePetsc, PETSC_NULL, PETSC_NULL);
    SlepcInitialize(&argPetscSize, &argCStylePetsc, PETSC_NULL, PETSC_NULL);

    delete[] argCStylePetsc;

    // Stop PETSc when error
    // PetscPushErrorHandler(PetscAbortErrorHandler, NULL);
    */
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

    // Delete Options
    delete option;
    option = NULL;

    // Finalize MPI, PETSc and SLEPc
    //PetscFinalize();
    //SlepcFinalize();
    //MPI_Finalize();

    finaOne = true;
  }
}

Options& SmallFem::getOptions(void){
  return *option;
}
