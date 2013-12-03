#include <fstream>

#include "FormulationTyped.h"
#include "Exception.h"
#include "System.h"

#include "SystemAbstract.h"


using namespace std;

SystemAbstract::~SystemAbstract(void){
}

void SystemAbstract::constraint(const Formulation& formulation){
  /*
  // Check the formulation type
  string formulationType = formulation.getType();
  string systemType      = getType();

  if(formulationType.compare(systemType) != 0)
    throw Exception("SystemAbstract::constraint -- scalar type mismatch");
  */
  // Check domaine
  const FunctionSpace& cFS = formulation.fs();
  if(&(cFS.getSupport().getMesh()) != &(fs->getSupport().getMesh()))
    throw
      Exception
      ("SystemAbstract::constraint -- FunctionSpaces must use the same mesh");

  // Scalar of Vector ?
  if(cFS.isScalar() != fs->isScalar())
    throw Exception("SystemAbstract::constraint -- Scalar/Vector mismatch");

  // Assemble and Solve the problem
  const FormulationTyped<double>& myFormulation =
    static_cast<const FormulationTyped<double>&>(formulation);

  System cSys(myFormulation);
  cSys.assemble();
  cSys.solve();

  // Fix This SystemAbstract Dofs with sys Solution //
  const vector<Dof> dof = cFS.getAllDofs();
  const size_t     nDof = dof.size();

  const DofManager&  cDofM = cSys.getDofManager();
  fullVector<double> cSol;
  cSys.getSolution(cSol);

  for(size_t i = 0; i < nDof; i++)
    dofM->fixValue(dof[i], cSol(cDofM.getGlobalId(dof[i])));
}

void SystemAbstract::writeMatrix(std::string fileName,
                                 std::string matrixName) const{
  ofstream stream;
  stream.open(fileName.c_str());
  stream << "writeMatrix not implemented for this system" << endl;
  stream.close();
}
