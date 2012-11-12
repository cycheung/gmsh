#include "EigenSystem.h"

using namespace std;

EigenSystem::EigenSystem(const EigenFormulation& eFormulation){
  // Get EigenFormulation //
  this->eFormulation = &eFormulation;
  this->fs           = &(eFormulation.fs());

  // Get Dof Manager //
  dofM = new DofManager();
  dofM->addDof(*fs);

  // Get DofManager Data //
  size = fs->dofNumber();
  
  // Is the Problem a General EigenValue Problem ? //
  isGeneral = eFormulation.isGeneral();

  // Create EigenSystem //
  // Linear System A
  linSysA = new linearSystemPETSc<double>();
  linSysA->allocate(size);

  // Linear System B
  if(isGeneral){
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

  // The EigenSystem is not assembled //
  isAssembled = false;
}

EigenSystem::~EigenSystem(void){
  if(eigenVector)
    delete eigenVector;
  
  if(eigenValue);
  delete eigenValue;

  if(eSys)
    delete eSys;

  delete linSysA;

  if(isGeneral)
    delete linSysB;

  delete dofM;
  // EigenSystem is not responsible for deleting 'Formulations'
}

void EigenSystem::assemble(void){
  // Get GroupOfDofs //
  const std::vector<GroupOfDof*>& group = fs->getAllGroups();
  const int E = fs->groupNumber();

  // Get Sparcity Pattern & PreAllocate//
  if(isGeneral)
    for(int i = 0; i < E; i++)
      sparcityGeneral(*(group[i]));  

  else
    for(int i = 0; i < E; i++)
      sparcity(*(group[i]));  

  linSysA->preAllocateEntries();
  
  if(isGeneral)
    linSysB->preAllocateEntries();

  // Assemble EigenSystem //
  if(isGeneral)
    for(int i = 0; i < E; i++)
      assembleGeneral(*(group[i]));  

  else
    for(int i = 0; i < E; i++)
      assemble(*(group[i]));  

  // The EigenSystem is assembled //
  isAssembled = true;  
}

void EigenSystem::fixDof(const GroupOfElement& goe, double value){
  const vector<const MElement*>&  element = goe.getAll();
  unsigned int                   nElement = goe.getNumber();
  
  for(unsigned int i = 0; i < nElement; i++){
    vector<Dof>         dof = fs->getKeys(*element[i]);
    const unsigned int nDof = dof.size();
    
    for(unsigned int j = 0; j < nDof; j++)
      dofM->fixValue(dof[j], value);
  }
}

void EigenSystem::solve(unsigned int nEigenValues){
  // Is the EigenSystem assembled ? //
  if(!isAssembled)
    assemble();

  // Solve //
  eSys = new EigenSolver(linSysA, linSysB);
  eSys->solve(nEigenValues);

  // Get Solution //
  nEigenValue = eSys->getNumEigenValues();
  
  eigenValue  = new vector<complex<double> >(nEigenValue);
  eigenVector = new vector<vector<complex<double> > >(nEigenValue);
  
  for(unsigned int i = 0; i < nEigenValue; i++)
    (*eigenValue)[i] = eSys->getEigenValue(i);

  for(unsigned int i = 0; i < nEigenValue; i++)
    (*eigenVector)[i] = eSys->getEigenVector(i);  
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
      linSysA->addToRightHandSide(dofI, fixed.second); 
    }
       
    else{
      // If unknown Dof
      for(int j = 0; j < N; j++){
	int dofJ = dofM->getGlobalId(*(dof[j]));

	linSysA->addToMatrix(dofI, dofJ, 
			     eFormulation->weakA(i, j, group));
      }
      
      linSysA->addToRightHandSide(dofI, 
				  eFormulation->rhs(i, group)); 
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
      linSysA->addToRightHandSide(dofI, fixed.second); 

      linSysB->addToMatrix(dofI, dofI, 1);
      linSysB->addToRightHandSide(dofI, fixed.second); 
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
      
      linSysA->addToRightHandSide(dofI, 
				  eFormulation->rhs(i, group)); 

      linSysB->addToRightHandSide(dofI, 
				  eFormulation->rhs(i, group)); 
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

