#include "BasisGenerator.h"
#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"

#include "FEMSolution.h"
#include "fullMatrix.h"
#include "DofManager.h"
#include "Mesh.h"

#include "Exception.h"
#include "SmallFem.h"

using namespace std;

void compute(const Options& option){
  // Get Domains //
  Mesh msh(option.getValue("-msh")[0]);
  GroupOfElement domain = msh.getFromPhysical(7);

  // Get FunctionSpace //
  const size_t   order = atoi(option.getValue("-o")[0].c_str());
  Basis*         basis;
  FunctionSpace* fSpace;

  if(option.getValue("-type")[0].compare("lagrange") == 0){
    // If Lagrange
    basis =
      BasisGenerator::generate(domain.get(0).getType(),
                               0, order, "lagrange");

    fSpace = new FunctionSpaceScalar(domain, *basis);
  }

  else if(option.getValue("-type")[0].compare("scalar") == 0){
    // If Scalar
    basis =
      BasisGenerator::generate(domain.get(0).getType(),
                               0, order, "hierarchical");

    fSpace = new FunctionSpaceScalar(domain, *basis);
  }

  else if(option.getValue("-type")[0].compare("vector") == 0){
    // If Vector
    basis =
      BasisGenerator::generate(domain.get(0).getType(),
                               1, order, "hierarchical");

    fSpace = new FunctionSpaceVector(domain, *basis);
  }

  else
    throw Exception("Unknown FunctionSpace type: %s",
                    option.getValue("-type")[0].c_str());

  // Enumerate Dofs //
  DofManager dofM;
  dofM.addToDofManager(fSpace->getAllGroups());
  dofM.generateGlobalIdSpace();

  // FunctionSpace is solution of FEM problem                   //
  // For every global basis, a vector with all zero excepte one //
  // (associated to the basis) can be created                   //

  // Do this with FEMSolution //
  const size_t nCoef = fSpace->dofNumber();
  fullVector<double> coef(nCoef);
  FEMSolution sol;

  // Init coefs to 0
  for(size_t i = 0; i < nCoef; i++)
    coef(i) = 0;

  // Set each coef to one (one at a time)
  for(size_t i = 0; i < nCoef; i++){
    coef(i) = 1;
    sol.addCoefficients(i, i, *fSpace, dofM, coef);
    coef(i) = 0;
  }

  // Write //
  sol.write("function_space");

  // Clean //
  delete fSpace;
  delete basis;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-type");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
