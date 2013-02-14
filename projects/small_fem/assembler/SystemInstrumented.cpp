#include "FormulationProjectionScalar.h"
#include "FormulationProjectionVector.h"
#include "BasisGenerator.h"
#include "BasisLocal.h"
#include "Exception.h"

#include "Timer.h"
#include "SystemInstrumented.h"

using namespace std;

SystemInstrumented::SystemInstrumented(const Formulation& formulation){
  // Get Formulation //
  this->formulation = &formulation;
  this->fs          = &(formulation.fs());

  // Get Dof Manager //
  dofM = new DofManager();
  dofM->addToGlobalIdSpace(fs->getAllGroups());

  // Get DofManager Data //
  size = fs->dofNumber();

  // Create System //
  x      = new fullVector<double>(size);
  linSys = new linearSystemPETSc<double>();
  linSys->allocate(size);

  // The system is not assembled and not solved //
  assembled = false;
  solved    = false;

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
  delete x;
  delete linSys;
  delete dofM;
  // System is not responsible for deleting 'Formulations'
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

void SystemInstrumented::fixCoef(const GroupOfElement& goe, double value){
  const vector<pair<const MElement*, ElementData> >&
    element = goe.getAll();

  const unsigned int nElement = goe.getNumber();

  for(unsigned int i = 0; i < nElement; i++){
    vector<Dof>         dof = fs->getKeys(*element[i].first);
    const unsigned int nDof = dof.size();

    for(unsigned int j = 0; j < nDof; j++)
      dofM->fixValue(dof[j], value);
  }
}

void SystemInstrumented::dirichlet(GroupOfElement& goe,
		       double (*f)(fullVector<double>& xyz)){

  // Check if Scalar Problem //
  if(!fs->isScalar())
    throw Exception("Cannot impose Vectorial Dirichlet Conditions on a Scalar Problem");

  // New FunctionSpace, on the Dirichlet Domain: dirFS //
  // WARNING: The support of the dirFS *MUST* have the fs Mesh
  //  --> So we have the same Dof Numbering

  if(&(goe.getMesh()) != &(fs->getSupport().getMesh()))
    throw Exception("Dirichlet Domain must come from the FunctionSpace Domain's Mesh");

  BasisLocal* dirBasis = BasisGenerator::generate(goe.get(0).getType(),
                                                  fs->getBasis(0).getType(),
                                                  fs->getBasis(0).getOrder(),
                                                  "hierarchical");
  FunctionSpaceScalar dirFS(goe, *dirBasis);

  // Solve The Projection Of f on the Dirichlet Domain with dirFS //
  FormulationProjectionScalar projection(f, dirFS);
  SystemInstrumented sysProj(projection);

  sysProj.assemble();
  sysProj.solve();

  // Fix This System Dofs with sysProj Solution //
  const vector<const Dof*> dof = dirFS.getAllDofs();
  const unsigned int      nDof = dof.size();

  const DofManager&        dirDofM = sysProj.getDofManager();
  const fullVector<double>& dirSol = sysProj.getSol();

  for(unsigned int i = 0; i < nDof; i++)
    dofM->fixValue(*dof[i], dirSol(dirDofM.getGlobalId(*dof[i])));

  delete dirBasis;
}

void SystemInstrumented::dirichlet(GroupOfElement& goe,
		       fullVector<double> (*f)(fullVector<double>& xyz)){

  // Check if Scalar Problem //
  if(fs->isScalar())
    throw Exception("Cannot impose Scalar Dirichlet Conditions on a Vectorial Problem");

  // New FunctionSpace, on the Dirichlet Domain: dirFS //
  // WARNING: The support of the dirFS *MUST* have the fs Mesh
  //  --> So we have the same Dof Numbering

  if(&(goe.getMesh()) != &(fs->getSupport().getMesh()))
    throw Exception("Dirichlet Domain must come from the FunctionSpace Domain's Mesh");

  BasisLocal* dirBasis = BasisGenerator::generate(goe.get(0).getType(),
                                                  fs->getBasis(0).getType(),
                                                  fs->getBasis(0).getOrder(),
                                                  "hierarchical");

  FunctionSpaceVector dirFS(goe, *dirBasis);

  // Solve The Projection Of f on the Dirichlet Domain with dirFS //
  FormulationProjectionVector projection(f, dirFS);
  SystemInstrumented sysProj(projection);

  sysProj.assemble();
  sysProj.solve();

  // Fix This System Dofs with sysProj Solution //
  const vector<const Dof*> dof = dirFS.getAllDofs();
  const unsigned int      nDof = dof.size();

  const DofManager&        dirDofM = sysProj.getDofManager();
  const fullVector<double>& dirSol = sysProj.getSol();

  for(unsigned int i = 0; i < nDof; i++)
    dofM->fixValue(*dof[i], dirSol(dirDofM.getGlobalId(*dof[i])));

  delete dirBasis;
}

void SystemInstrumented::solve(void){
  // Is the System assembled ? //
  if(!assembled)
    assemble();

  // Solve //
  linSys->systemSolve();

  // Write Sol
  double xi;

  for(int i = 0; i < size; i++){
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

