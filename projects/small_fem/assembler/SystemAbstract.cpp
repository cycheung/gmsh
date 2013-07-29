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
  const vector<const MElement*>& element = goe.getAll();
  const size_t                  nElement = goe.getNumber();

  for(size_t i = 0; i < nElement; i++){
    vector<Dof>   dof = fs->getKeys(*element[i]);
    const size_t nDof = dof.size();

    for(size_t j = 0; j < nDof; j++)
      dofM->fixValue(dof[j], value);
  }
}

void SystemAbstract::
dirichlet(GroupOfElement& goe,
          double (*f)(fullVector<double>& xyz)){

  // Check if Scalar Problem //
  if(!fs->isScalar())
    throw
      Exception
      ("Cannot impose Vectorial Dirichlet Conditions on a Scalar Problem");

  // New FunctionSpace, on the Dirichlet Domain: dirFS //
  // WARNING: The support of the dirFS *MUST* have the fs Mesh
  //  --> So we have the same Dof Numbering

  if(&(goe.getMesh()) != &(fs->getSupport().getMesh()))
    throw
      Exception
      ("Dirichlet Domain must come from the FunctionSpace Domain's Mesh");

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
  const vector<Dof> dof = dirFS.getAllDofs();
  const size_t     nDof = dof.size();

  const DofManager&        dirDofM = sysProj.getDofManager();
  const fullVector<double>& dirSol = sysProj.getSol();

  for(size_t i = 0; i < nDof; i++)
    dofM->fixValue(dof[i], dirSol(dirDofM.getGlobalId(dof[i])));

  delete dirBasis;
}

void SystemAbstract::
dirichlet(GroupOfElement& goe,
          fullVector<double> (*f)(fullVector<double>& xyz)){

  // Check if Scalar Problem //
  if(fs->isScalar())
    throw
      Exception
      ("Cannot impose Scalar Dirichlet Conditions on a Vectorial Problem");

  // New FunctionSpace, on the Dirichlet Domain: dirFS //
  // WARNING: The support of the dirFS *MUST* have the fs Mesh
  //  --> So we have the same Dof Numbering

  if(&(goe.getMesh()) != &(fs->getSupport().getMesh()))
    throw
      Exception
      ("Dirichlet Domain must come from the FunctionSpace Domain's Mesh");

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
  const vector<Dof> dof = dirFS.getAllDofs();
  const size_t     nDof = dof.size();

  const DofManager&        dirDofM = sysProj.getDofManager();
  const fullVector<double>& dirSol = sysProj.getSol();

  for(size_t i = 0; i < nDof; i++)
    dofM->fixValue(dof[i], dirSol(dirDofM.getGlobalId(dof[i])));

  delete dirBasis;
}

void SystemAbstract::assemble(Mat& A,
                              Vec& b,
                              size_t elementId,
                              const GroupOfDof& group,
                              formulationPtr& term){

  const vector<Dof>& dof = group.getDof();
  const size_t N = group.size();

  size_t dofI;
  size_t dofJ;

  //PetscInt petscI;
  //PetscInt petscJ;
  //PetscScalar petscV;

  for(size_t i = 0; i < N; i++){
    dofI = dofM->getGlobalId(dof[i]);

    // If not a fixed Dof line: assemble
    if(dofI != DofManager::isFixedId()){
      for(size_t j = 0; j < N; j++){
        dofJ = dofM->getGlobalId(dof[j]);

        // If not a fixed Dof
        if(dofJ != DofManager::isFixedId()){
          MatSetValue(A, dofI, dofJ,
                      (formulation->*term)(i, j, elementId), ADD_VALUES);

          // TODO: Concider using MatSetValueS
        }

        // If fixed Dof (for column 'dofJ'):
        //    add to right hand side (with a minus sign) !
        else{
            VecSetValue(b, dofI,
                        -1 * dofM->getValue(dof[j]) *
                        (formulation->*term)(i, j, elementId),
                        ADD_VALUES);

            // TODO: Concider using VecSetValueS
        }
      }

      VecSetValue(b, dofI,
                  formulation->rhs(i, elementId),
                  ADD_VALUES);

      // TODO: Concider using VecSetValueS
    }
  }
}

void SystemAbstract::assemble(Mat& A,
                              size_t elementId,
                              const GroupOfDof& group,
                              formulationPtr& term){

  const vector<Dof>& dof = group.getDof();
  const size_t N = group.size();

  size_t dofI;
  size_t dofJ;

  //PetscInt petscI;
  //PetscInt petscJ;
  //PetscScalar petscV;

  for(size_t i = 0; i < N; i++){
    dofI = dofM->getGlobalId(dof[i]);

    // If not a fixed Dof line: assemble
    if(dofI != DofManager::isFixedId()){
      for(size_t j = 0; j < N; j++){
        dofJ = dofM->getGlobalId(dof[j]);

        // If not a fixed Dof
        if(dofJ != DofManager::isFixedId()){
          MatSetValue(A, dofI, dofJ,
                      (formulation->*term)(i, j, elementId), ADD_VALUES);

          // TODO: Concider using MatSetValueS
        }
      }
    }
  }
}

void SystemAbstract::sparsity(PetscInt* nonZero,
                              UniqueSparsity& uniqueSparsity,
                              const GroupOfDof& group){

  const vector<Dof>& dof = group.getDof();
  const size_t N = group.size();

  size_t dofI;
  size_t dofJ;

  // Add each column only one

  for(size_t i = 0; i < N; i++){
    dofI = dofM->getGlobalId(dof[i]);

    // Add non fixed Dof
    if(dofI != DofManager::isFixedId()){
      for(size_t j = 0; j < N; j++){
        dofJ = dofM->getGlobalId(dof[j]);

        // Add non fixed Dof
        if(dofJ != DofManager::isFixedId()){
          // Check if pair (dofI, dofJ) allready inserted
          if(uniqueSparsity.empty() ||
             uniqueSparsity[dofI].insert(dofJ).second)
              nonZero[dofI]++;
        }
      }
    }
  }
}
