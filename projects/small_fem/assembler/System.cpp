#include "FormulationProjectionScalar.h"
#include "FormulationProjectionVector.h"
#include "Exception.h"

#include "System.h"

using namespace std;

System::System(const Formulation& formulation){
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
}

System::~System(void){
  delete x;
  delete linSys;
  delete dofM;
  // System is not responsible for deleting 'Formulations'
}

void System::assemble(void){
  // Get GroupOfDofs //
  const vector<GroupOfDof*>& group = fs->getAllGroups();
  const int E = fs->groupNumber();

  // Set to put Fixed Dof only ones 
  // (cannot use both  setValue and add Value
  //  in PETSc)
  fixedOnes = new set<const Dof*, DofComparator>();

  // Get Sparcity Pattern & PreAllocate//
  for(int i = 0; i < E; i++)
    sparcity(*(group[i]));  

  linSys->preAllocateEntries();

  // Assemble System //
  for(int i = 0; i < E; i++)
    assemble(*(group[i]));  

  // The system is assembled //
  delete fixedOnes;
  assembled = true;  
}

void System::fixCoef(const GroupOfElement& goe, double value){
  const vector<const MElement*>&  element = goe.getAll();
  unsigned int                   nElement = goe.getNumber();
  
  for(unsigned int i = 0; i < nElement; i++){
    vector<Dof>         dof = fs->getKeys(*element[i]);
    const unsigned int nDof = dof.size();
    
    for(unsigned int j = 0; j < nDof; j++)
      dofM->fixValue(dof[j], value);
  }
}

void System::dirichlet(const GroupOfElement& goe, 
		       double (*f)(fullVector<double>& xyz)){

  // Check if Scalar Problem //
  if(!fs->isScalar())
    throw Exception("Cannot impose Vectorial Dirichlet Conditions on a Scalar Problem");

  // New FunctionSpace, on the Dirichlet Domain: dirFS //
  // WARNING: The support of the dirFS *MUST* have the fs Mesh
  //  --> So we have the same Dof Numbering

  if(&(goe.getMesh()) != &(fs->getSupport().getMesh()))
    throw Exception("Dirichlet Domain must come from the FunctionSpace Domain's Mesh");

  FunctionSpaceNode dirFS(goe, fs->getOrder());

  // Solve The Projection Of f on the Dirichlet Domain with dirFS //
  FormulationProjectionScalar projection(f, dirFS);
  System sysProj(projection);

  sysProj.assemble();
  sysProj.solve();

  // Fix This System Dofs with sysProj Solution //
  const vector<const Dof*> dof = dirFS.getAllDofs();
  const unsigned int      nDof = dof.size();

  const DofManager&        dirDofM = sysProj.getDofManager();
  const fullVector<double>& dirSol = sysProj.getSol();

  for(unsigned int i = 0; i < nDof; i++)
    dofM->fixValue(*dof[i], dirSol(dirDofM.getGlobalId(*dof[i]))); 
}

void System::dirichlet(const GroupOfElement& goe, 
		       fullVector<double> (*f)(fullVector<double>& xyz)){

  // Check if Scalar Problem //
  if(fs->isScalar())
    throw Exception("Cannot impose Scalar Dirichlet Conditions on a Vectorial Problem");

  // New FunctionSpace, on the Dirichlet Domain: dirFS //
  // WARNING: The support of the dirFS *MUST* have the fs Mesh
  //  --> So we have the same Dof Numbering

  if(&(goe.getMesh()) != &(fs->getSupport().getMesh()))
    throw Exception("Dirichlet Domain must come from the FunctionSpace Domain's Mesh");

  FunctionSpaceEdge dirFS(goe, fs->getOrder());

  // Solve The Projection Of f on the Dirichlet Domain with dirFS //
  FormulationProjectionVector projection(f, dirFS);
  System sysProj(projection);

  sysProj.assemble();
  sysProj.solve();

  // Fix This System Dofs with sysProj Solution //
  const vector<const Dof*> dof = dirFS.getAllDofs();
  const unsigned int      nDof = dof.size();

  const DofManager&        dirDofM = sysProj.getDofManager();
  const fullVector<double>& dirSol = sysProj.getSol();

  for(unsigned int i = 0; i < nDof; i++)
    dofM->fixValue(*dof[i], dirSol(dirDofM.getGlobalId(*dof[i]))); 
}


void System::solve(void){
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

void System::assemble(GroupOfDof& group){
  const vector<const Dof*>& dof = group.getAll();
  const int N = group.getNumber();

  for(int i = 0; i < N; i++){
    pair<bool, double> fixed = dofM->getValue(*(dof[i]));
    int dofI = dofM->getGlobalId(*(dof[i]));

    if(fixed.first){
      // If fixed Dof
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
	int dofJ = dofM->getGlobalId(*(dof[j]));

	linSys->addToMatrix(dofI, dofJ, 
			    formulation->weak(i, j, group));
      }
      
      linSys->addToRightHandSide(dofI, 
				 formulation->rhs(i, group)); 
    }
  } 
}

void System::sparcity(GroupOfDof& group){
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

