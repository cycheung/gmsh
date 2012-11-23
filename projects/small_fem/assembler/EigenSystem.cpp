#include "Exception.h"
#include "System.h"
#include "FormulationProjectionScalar.h"
#include "EigenSystem.h"

using namespace std;

EigenSystem::EigenSystem(const EigenFormulation& eFormulation){
  // Get EigenFormulation //
  this->eFormulation = &eFormulation;
  this->fs           = &(eFormulation.fs());

  // Get Dof Manager //
  dofM = new DofManager();
  dofM->addToGlobalIdSpace(fs->getAllGroups());

  // Get DofManager Data //
  size = fs->dofNumber();
  
  // Is the Problem a General EigenValue Problem ? //
  general = eFormulation.isGeneral();

  // Create EigenSystem //
  // Linear System A
  linSysA = new linearSystemPETSc<double>();
  linSysA->allocate(size);

  // Linear System B
  if(general){
    linSysB = new linearSystemPETSc<double>();
    linSysB->allocate(size);
  }
  
  else{
    linSysB = NULL;  
  }

  // eSys will be created at solving point
  eSys        = NULL; 
  eigenValue  = NULL;
  eigenVector = NULL;

  // The EigenSystem is not assembled and not solved//
  nEigenValue = 0;
  assembled   = false;
  solved      = false;
}

EigenSystem::~EigenSystem(void){
  if(eigenVector)
    delete eigenVector;
  
  if(eigenValue)
    delete eigenValue;

  if(eSys)
    delete eSys;

  delete linSysA;

  if(general)
    delete linSysB;

  delete dofM;
  // EigenSystem is not responsible for deleting 'Formulations'
}

void EigenSystem::assemble(void){
  // Get GroupOfDofs //
  const vector<GroupOfDof*>& group = fs->getAllGroups();
  const int E = fs->groupNumber();

  // Get Sparcity Pattern & PreAllocate//
  if(general)
    for(int i = 0; i < E; i++)
      sparcityGeneral(*(group[i]));  

  else
    for(int i = 0; i < E; i++)
      sparcity(*(group[i]));  

  linSysA->preAllocateEntries();
  
  if(general)
    linSysB->preAllocateEntries();

  // Assemble EigenSystem //
  if(general)
    for(int i = 0; i < E; i++)
      assembleGeneral(*(group[i]));  

  else
    for(int i = 0; i < E; i++)
      assemble(*(group[i]));  

  // The EigenSystem is assembled //
  assembled = true;  
}

void EigenSystem::fixCoef(const GroupOfElement& goe, double value){
  const vector<const MElement*>&  element = goe.getAll();
  unsigned int                   nElement = goe.getNumber();
  
  for(unsigned int i = 0; i < nElement; i++){
    vector<Dof>         dof = fs->getKeys(*element[i]);
    const unsigned int nDof = dof.size();
    
    for(unsigned int j = 0; j < nDof; j++)
      dofM->fixValue(dof[j], value);
  }
}

void EigenSystem::dirichlet(const GroupOfElement& goe, 
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

void EigenSystem::solve(unsigned int nEigenValues){
  // Check nEigenValues
  if(nEigenValues > size)
    throw 
      Exception("I cannot compute more Eigenvalues (%d) than the number of unknowns (%d)",
		nEigenValues, size);

  // Is the EigenSystem assembled ? //
  if(!assembled)
    assemble();

  // Solve //
  eSys = new EigenSolver(linSysA, linSysB, false);
  eSys->solve(nEigenValues, "smallest");
  
  // Get Solution //
  nEigenValue = eSys->getNumEigenValues();
  
  eigenValue  = new vector<complex<double> >(nEigenValue);
  eigenVector = new vector<vector<complex<double> > >(nEigenValue);
  
  for(unsigned int i = 0; i < nEigenValue; i++)
    (*eigenValue)[i] = eSys->getEigenValue(i);

  for(unsigned int i = 0; i < nEigenValue; i++)
    (*eigenVector)[i] = eSys->getEigenVector(i);  

  // System solved ! //
  solved = true;
}

void EigenSystem::assemble(GroupOfDof& group){
  const vector<const Dof*>& dof = group.getAll();
  const int N = group.getNumber();

  for(int i = 0; i < N; i++){
    pair<bool, double> fixed = dofM->getValue(*(dof[i]));
    int dofI = dofM->getGlobalId(*(dof[i]));

    if(fixed.first){
      // If fixed Dof
      linSysA->addToMatrix(dofI, dofI, 1);
    }
       
    else{
      // If unknown Dof
      for(int j = 0; j < N; j++){
	int dofJ = dofM->getGlobalId(*(dof[j]));

	linSysA->addToMatrix(dofI, dofJ, 
			     eFormulation->weakA(i, j, group));
      }
    }
  }
}

void EigenSystem::assembleGeneral(GroupOfDof& group){
  const vector<const Dof*>& dof = group.getAll();
  const int N = group.getNumber();

  for(int i = 0; i < N; i++){
    pair<bool, double> fixed = dofM->getValue(*(dof[i]));
    int dofI = dofM->getGlobalId(*(dof[i]));

    if(fixed.first){
      // If fixed Dof
      linSysA->addToMatrix(dofI, dofI, 1);
      linSysB->addToMatrix(dofI, dofI, 1);
    }
       
    else{
      // If unknown Dof
      for(int j = 0; j < N; j++){
	int dofJ = dofM->getGlobalId(*(dof[j]));

	linSysA->addToMatrix(dofI, dofJ, 
			     eFormulation->weakA(i, j, group));

	linSysB->addToMatrix(dofI, dofJ, 
			     eFormulation->weakB(i, j, group));
      }
    }
  } 
}

void EigenSystem::sparcity(GroupOfDof& group){
  const vector<const Dof*>& dof = group.getAll();
  const int N = group.getNumber();

  for(int i = 0; i < N; i++){
    pair<bool, double> fixed = dofM->getValue(*(dof[i]));
    int dofI = dofM->getGlobalId(*(dof[i]));

    if(fixed.first)
      // If fixed Dof
      linSysA->insertInSparsityPattern(dofI, dofI);
    
    else
      // If unknown Dof
      for(int j = 0; j < N; j++){
	int dofJ = dofM->getGlobalId(*(dof[j]));
	
	linSysA->insertInSparsityPattern(dofI, dofJ);
      } 
  } 
}

void EigenSystem::sparcityGeneral(GroupOfDof& group){
  const vector<const Dof*>& dof = group.getAll();
  const int N = group.getNumber();

  for(int i = 0; i < N; i++){
    pair<bool, double> fixed = dofM->getValue(*(dof[i]));
    int dofI = dofM->getGlobalId(*(dof[i]));

    if(fixed.first){
      // If fixed Dof
      linSysA->insertInSparsityPattern(dofI, dofI);
      linSysB->insertInSparsityPattern(dofI, dofI);
    }

    else{
      // If unknown Dof
      for(int j = 0; j < N; j++){
	int dofJ = dofM->getGlobalId(*(dof[j]));

	linSysA->insertInSparsityPattern(dofI, dofJ);
	linSysB->insertInSparsityPattern(dofI, dofJ);
      } 
    }
  } 
}

