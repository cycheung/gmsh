#include "Timer.h"
#include "DofFixedException.h"
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
  // Enumerate //
  dofM->generateGlobalIdSpace();

  // Init System //
  x      = new fullVector<double>(dofM->getDofNumber());
  linSys = new linearSystemPETSc<double>();
  linSys->allocate(dofM->getDofNumber());

  // Get Groups //
  const unsigned int E = fs->getSupport().getNumber();
  const vector<GroupOfDof*>& group = fs->getAllGroups();

  // Get Sparsity Pattern & PreAllocate//
  Timer timer;
  timer.start();

  for(unsigned int i = 0; i < E; i++)
    SystemAbstract::sparsity(*linSys, *group[i]);

  linSys->preAllocateEntries();

  timer.stop();
  preAlloc = timer.time();

  // Assemble System //
  for(unsigned int i = 0; i < E; i++)
    assemble(*group[i], i);

  // The system is assembled //
  assembled = true;
}

void SystemInstrumented::assemble(const GroupOfDof& group,
                                  unsigned int elementId){
  Timer timer;
  double a;
  double b;

  unsigned int dofI;
  unsigned int dofJ;

  const vector<const Dof*>& dof = group.getAll();
  const unsigned int N = group.getNumber();

  for(unsigned int i = 0; i < N; i++){
    try{
      timer.start();
      dofI = dofM->getGlobalId(*(dof[i]));
      timer.stop();

      dofLookTime += timer.time();
      dofLookCall++;

      for(unsigned int j = 0; j < N; j++){
        try{
          timer.start();
          dofJ = dofM->getGlobalId(*(dof[j]));
          timer.stop();

          dofLookTime += timer.time();
          dofLookCall++;

          timer.start();
          a = formulation->weak(i, j, elementId);
          timer.stop();

          totLHSTime += timer.time();
          totLHSCall++;

          timer.start();
          linSys->addToMatrix(dofI, dofJ, a);
          timer.stop();

          totAddLHSTime += timer.time();
          totAddLHSCall++;
        }

        catch(DofFixedException& fixedDof){
          // If fixed Dof (for column 'dofJ'):
          //    add to right hand side (with a minus sign) !
          timer.start();
          b = -1 * fixedDof.getValue() * (formulation->weak)(i, j, elementId);
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

      timer.start();
      b = formulation->rhs(i, elementId);
      timer.stop();

      totRHSTime += timer.time();
      totRHSCall++;

      timer.start();
      linSys->addToRightHandSide(dofI, b);
      timer.stop();

      totAddRHSTime += timer.time();
      totAddRHSCall++;
    }

    catch(DofFixedException& any){
      // Don't add fixed Dof (for line 'dofI')
    }
  }
}

