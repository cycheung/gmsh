#include "Timer.h"
#include "SystemInstrumented.h"

using namespace std;

SystemInstrumented::SystemInstrumented(const Formulation& formulation):
  System(formulation){

  // Timers //
  totLHSTime = 0;
  totLHSCall = 0;

  totAddLHSTime = 0;
  totAddLHSCall = 0;

  totRHSTime = 0;
  totRHSCall = 0;

  totAddRHSTime = 0;
  totAddRHSCall = 0;

  dofLookTime = 0;
  dofLookCall = 0;
}

SystemInstrumented::~SystemInstrumented(void){
}

void SystemInstrumented::assemble(void){
  // Get GroupOfDofs //
  const vector<GroupOfDof*>& group = fs->getAllGroups();
  const int E = fs->groupNumber();

  // Set to put Fixed Dof only ones
  // (cannot use both  setValue and add Value
  //  in PETSc)
  fixedOnes = new set<const Dof*, DofComparator>();

  // Get Sparsity Pattern & PreAllocate//
  Timer timer;
  timer.start();

  for(int i = 0; i < E; i++)
    sparsity(*(group[i]));

  linSys->preAllocateEntries();

  timer.stop();
  preAlloc = timer.time();

  // Assemble System //
  for(int i = 0; i < E; i++)
    assemble(*(group[i]));

  // The system is assembled //
  delete fixedOnes;
  assembled = true;
}

void SystemInstrumented::solve(void){
  // Is the System assembled ? //
  if(!assembled)
    assemble();

  // Solve //
  linSys->systemSolve();

  // Write Sol
  const unsigned int size = fs->dofNumber();
  double xi;

  for(unsigned int i = 0; i < size; i++){
    linSys->getFromSolution(i, xi);
    (*x)(i) = xi;
  }

  // System solved ! //
  solved = true;
}

void SystemInstrumented::assemble(GroupOfDof& group){
  Timer timer;
  double a;
  double b;

  const vector<const Dof*>& dof = group.getAll();
  const int N = group.getNumber();

  for(int i = 0; i < N; i++){
    timer.start();
    pair<bool, double> fixed = dofM->getValue(*(dof[i]));
    int dofI = dofM->getGlobalId(*(dof[i]));
    timer.stop();

    dofLookTime += timer.time();
    dofLookCall++;

    if(fixed.first){
      pair<
	set<const Dof*, DofComparator>::iterator,
	bool> ones = fixedOnes->insert(dof[i]);

      if(ones.second){
	linSys->addToMatrix(dofI, dofI, 1);
	linSys->addToRightHandSide(dofI, fixed.second);
      }
   }

    else{
      // If unknown Dof
      for(int j = 0; j < N; j++){
        timer.start();
	int dofJ = dofM->getGlobalId(*(dof[j]));
        timer.stop();

        dofLookTime += timer.time();
        dofLookCall++;

        timer.start();
        a = formulation->weak(i, j, group);
        timer.stop();

        totLHSTime += timer.time();
        totLHSCall++;

        timer.start();
	linSys->addToMatrix(dofI, dofJ, a);
        timer.stop();

        totAddLHSTime += timer.time();
        totAddLHSCall++;
      }

      timer.start();
      b = formulation->rhs(i, group);
      timer.stop();

      totRHSTime += timer.time();
      totRHSCall++;

      timer.start();
      linSys->addToRightHandSide(dofI, b);
      timer.stop();

      totAddRHSTime += timer.time();
      totAddRHSCall++;
    }
  }
}

void SystemInstrumented::sparsity(GroupOfDof& group){
  const vector<const Dof*>& dof = group.getAll();
  const int N = group.getNumber();

  for(int i = 0; i < N; i++){
    pair<bool, double> fixed = dofM->getValue(*(dof[i]));
    int dofI = dofM->getGlobalId(*(dof[i]));

    if(fixed.first){
      // If fixed Dof
      linSys->insertInSparsityPattern(dofI, dofI);
    }

    else{
      // If unknown Dof
      for(int j = 0; j < N; j++){
	int dofJ = dofM->getGlobalId(*(dof[j]));

	linSys->insertInSparsityPattern(dofI, dofJ);
      }
    }
  }
}

