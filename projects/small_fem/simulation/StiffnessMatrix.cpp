#include "Mesh.h"
#include "System.h"

#include "FormulationLaplace.h"
#include "FormulationCurl.h"

#include "Timer.h"
#include "SmallFem.h"

using namespace std;

void compute(const Options& option){
  // Get Domains //
  Mesh msh(option.getValue("-msh")[0]);
  GroupOfElement domain = msh.getFromPhysical(7);

  if(domain.getNumber() != 1)
    throw
      Exception
      ("To comute the Stiffness Matrix, I need only one element (%d given)",
       domain.getNumber());

  // Get Parameters //
  size_t order = atoi(option.getValue("-o")[0].c_str());

  // Assemble the asked matrix //
  string problemName = option.getValue("-type")[0];

  Formulation* formulation;

  if(!problemName.compare("grad"))
    formulation = new FormulationLaplace(domain, order);

  else if(!problemName.compare("curl"))
    formulation = new FormulationCurl(domain, order);

  else
    throw
      Exception
      ("I do not know what type of matrix '%s' is\nI know:%s%s\n",
       problemName.c_str(),
       " grad,"
       " curl ");

  System sys(*formulation);

  sys.assemble();

  // Write Matrix //
  string matrixName = option.getValue("-name")[0];
  string fileName   = matrixName;
  fileName.append(".m");

  sys.writeMatrix(fileName, matrixName);

  // Delete //
  delete formulation;
}

int main(int argc, char** argv){
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
