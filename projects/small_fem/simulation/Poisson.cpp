#include <iostream>
#include <sstream>

#include "Mesh.h"
#include "System.h"
#include "SystemHelper.h"

#include "FormulationPoisson.h"

#include "SmallFem.h"

using namespace std;

double fDirichlet0(fullVector<double>& xyz){
  return 0;
}

double fDirichlet1(fullVector<double>& xyz){
  return 0;
}

double fSource(fullVector<double>& xyz){
  return 1;
}

void compute(const Options& option){
  // Get Domains //
  Mesh msh(option.getValue("-msh")[0]);
  GroupOfElement    domain = msh.getFromPhysical(7);
  GroupOfElement boundary0 = msh.getFromPhysical(6);
  GroupOfElement boundary1 = msh.getFromPhysical(5);

  cout << "Number of Element: " << domain.getNumber()
       << endl << flush;

  // Get Order //
  size_t order = atoi(option.getValue("-o")[0].c_str());

  // Compute //
  FormulationPoisson poisson(domain, fSource, order);
  System sysPoisson(poisson);

  SystemHelper<double>::dirichlet(sysPoisson, boundary0, fDirichlet0);
  SystemHelper<double>::dirichlet(sysPoisson, boundary1, fDirichlet1);

  cout << "Poisson -- Order " << order
       << ": " << sysPoisson.getSize()
       << endl << flush;

  sysPoisson.assemble();
  sysPoisson.solve();

  // Write Sol //
  if(!option.getValue("-nopos").size()){
    FEMSolution feSol;
    sysPoisson.getSolution(feSol);
    feSol.write("poisson");
  }
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
