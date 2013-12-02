#include <fstream>

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

  const DofManager&  dirDofM = sysProj.getDofManager();
  fullVector<double> dirSol;
  sysProj.getSolution(dirSol);

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

  const DofManager&  dirDofM = sysProj.getDofManager();
  fullVector<double> dirSol;
  sysProj.getSolution(dirSol);

  for(size_t i = 0; i < nDof; i++)
    dofM->fixValue(dof[i], dirSol(dirDofM.getGlobalId(dof[i])));

  delete dirBasis;
}

void SystemAbstract::writeMatrix(std::string fileName,
                                 std::string matrixName) const{
  ofstream stream;
  stream.open(fileName.c_str());
  stream << "writeMatrix not implemented for this system" << endl;
  stream.close();
}
