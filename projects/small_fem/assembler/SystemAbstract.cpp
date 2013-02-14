#include "FormulationProjectionScalar.h"
#include "FormulationProjectionVector.h"
#include "BasisGenerator.h"
#include "BasisLocal.h"
#include "Exception.h"

#include "SystemAbstract.h"
#include "System.h"

using namespace std;

SystemAbstract::~SystemAbstract(void){
}

void SystemAbstract::fixCoef(const GroupOfElement& goe, double value){
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

void SystemAbstract::dirichlet(GroupOfElement& goe,
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

  delete dirBasis;
}

void SystemAbstract::dirichlet(GroupOfElement& goe,
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

  delete dirBasis;
}

void SystemAbstract::assemble(linearSystemPETSc<double>& sys,
                              const GroupOfDof& group,
                              formulationPtr& term){

  const vector<const Dof*>& dof = group.getAll();
  const unsigned int N = group.getNumber();

  pair<bool, double> fixed;
  unsigned int dofI;
  unsigned int dofJ;

  for(unsigned int i = 0; i < N; i++){
    fixed = dofM->getValue(*(dof[i]));
    dofM->getGlobalId(*(dof[i]));

    if(fixed.first){
      // If known Dof
      sys.addToMatrix(dofI, dofI, 1);
      sys.addToRightHandSide(dofI, fixed.second);
    }

    else{
      // If unknown Dof
      for(unsigned int j = 0; j < N; j++){
	dofM->getGlobalId(*(dof[j]));

	sys.addToMatrix(dofI, dofJ, (formulation->*term)(i, j, group));
      }

      sys.addToRightHandSide(dofI, formulation->rhs(i, group));
    }
  }
}

void SystemAbstract::sparsity(linearSystemPETSc<double>& sys,
                              const GroupOfDof& group){

  const vector<const Dof*>& dof = group.getAll();
  const unsigned int N = group.getNumber();

  pair<bool, double> fixed;
  unsigned int dofI;
  unsigned int dofJ;

  for(unsigned int i = 0; i < N; i++){
    fixed = dofM->getValue(*(dof[i]));
    dofI = dofM->getGlobalId(*(dof[i]));

    if(fixed.first){
      // If fixed Dof
      sys.insertInSparsityPattern(dofI, dofI);
    }

    else{
      // If unknown Dof
      for(unsigned int j = 0; j < N; j++){
	dofJ = dofM->getGlobalId(*(dof[j]));

	sys.insertInSparsityPattern(dofI, dofJ);
      }
    }
  }
}
