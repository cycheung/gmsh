#include "SystemShowFunctionSpace.h"

using namespace std;

SystemShowFunctionSpace::SystemShowFunctionSpace(const FunctionSpace& fs,
                                                 size_t functionNumber){
  // Get Formulation //
  this->formulation = NULL;
  this->fs          = &fs;

  // Get Dof Manager //
  dofM = new DofManager();
  dofM->addToDofManager(fs.getAllGroups());

  // Enumerate //
  dofM->generateGlobalIdSpace();

  // Alloc //
  A       = NULL;
  b       = NULL;
  //xPetsc  = NULL;
  x       = new fullVector<double>(dofM->getUnfixedDofNumber());
  fNumber = functionNumber;

  // The system is not assembled and not solved //
  assembled = false;
  solved    = false;
}

SystemShowFunctionSpace::~SystemShowFunctionSpace(void){
  // Done in System !
}

void SystemShowFunctionSpace::assemble(void){
  // The system is assembled //
  assembled = true;
}

void SystemShowFunctionSpace::solve(void){
  // Is the System assembled ? //
  if(!assembled)
    assemble();

  // Write Sol
  const size_t size = dofM->getUnfixedDofNumber();

  for(size_t i = 0; i < size; i++)
    if(i == fNumber)
      (*x)(i) = 1;

    else
      (*x)(i) = 0;

  // System solved ! //
  solved = true;
}
