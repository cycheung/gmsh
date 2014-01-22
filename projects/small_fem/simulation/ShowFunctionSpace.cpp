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
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement domain =
    msh.getFromPhysical(atoi(option.getValue("-phys")[1].c_str()));

  // Get FunctionSpace //
  const size_t   order = atoi(option.getValue("-o")[1].c_str());
  Basis*         basis;
  FunctionSpace* fSpace;

  if(option.getValue("-type")[1].compare("lagrange") == 0){
    // If Lagrange
    basis =
      BasisGenerator::generate(domain.get(0).getType(),
                               0, order, "lagrange");

    fSpace = new FunctionSpaceScalar(domain, *basis);
  }

  else if(option.getValue("-type")[1].compare("scalar") == 0){
    // If Scalar
    basis =
      BasisGenerator::generate(domain.get(0).getType(),
                               0, order, "hierarchical");

    fSpace = new FunctionSpaceScalar(domain, *basis);
  }

  else if(option.getValue("-type")[1].compare("vector") == 0){
    // If Vector
    basis =
      BasisGenerator::generate(domain.get(0).getType(),
                               1, order, "hierarchical");

    fSpace = new FunctionSpaceVector(domain, *basis);
  }

  else
    throw Exception("Unknown FunctionSpace type: %s",
                    option.getValue("-type")[1].c_str());

  // Enumerate Dofs //
  DofManager<double> dofM;
  dofM.addToDofManager(fSpace->getAllGroups());
  dofM.generateGlobalIdSpace();

  // FunctionSpace is solution of FEM problem                   //
  // For every global basis, a vector with all zero excepte one //
  // (associated to the basis) can be created                   //

  // Do this with FEMSolution //
  const size_t nCoef = fSpace->dofNumber();
  fullVector<double> coef(nCoef);
  FEMSolution<double> sol;

  // Init coefs to 0
  for(size_t i = 0; i < nCoef; i++)
    coef(i) = 0;

  // Set each coef to one (one at a time)
  size_t startCoef = 0;
  size_t stopCoef  = nCoef;

  try{
    if(option.getValue("-start").size() > 1)
      startCoef = atoi(option.getValue("-start")[1].c_str());
  }
  catch(...){
  }

  try{
    if(option.getValue("-stop").size() > 1)
      stopCoef = atoi(option.getValue("-stop")[1].c_str());
  }
  catch(...){
  }

  for(size_t i = startCoef; i < stopCoef; i++){
    cout << "# Dof: " << i << "/" << nCoef  << endl;
    coef(i) = 1;
    sol.addCoefficients(i - startCoef, i - startCoef, *fSpace, dofM, coef);
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
  SmallFem::Keywords("-msh,-o,-type,-start,-stop,-phys");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
