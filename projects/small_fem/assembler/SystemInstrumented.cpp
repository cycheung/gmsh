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
  const unsigned int E = fs->groupNumber();

  // Get Sparsity Pattern & PreAllocate//
  Timer timer;
  timer.start();

  for(unsigned int i = 0; i < E; i++)
    AbstractSystem::sparsity(*linSys, *(group[i]));

  linSys->preAllocateEntries();

  timer.stop();
  preAlloc = timer.time();

  // Assemble System //
  for(unsigned int i = 0; i < E; i++)
    assemble(*(group[i]));

  // The system is assembled //
  assembled = true;
}

void SystemInstrumented::assemble(GroupOfDof& group){
  Timer timer;
  double a;
  double b;

  const vector<const Dof*>& dof = group.getAll();
  const unsigned int N = group.getNumber();

  for(unsigned int i = 0; i < N; i++){
    timer.start();
    pair<bool, double> fixed = dofM->getValue(*(dof[i]));
    unsigned int dofI = dofM->getGlobalId(*(dof[i]));
    timer.stop();

    dofLookTime += timer.time();
    dofLookCall++;

    if(fixed.first){
      linSys->addToMatrix(dofI, dofI, 1);
      linSys->addToRightHandSide(dofI, fixed.second);
   }

    else{
      // If unknown Dof
      for(unsigned int j = 0; j < N; j++){
        timer.start();
	unsigned int dofJ = dofM->getGlobalId(*(dof[j]));
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

