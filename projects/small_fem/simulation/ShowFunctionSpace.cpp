#include <sstream>

#include "BasisGenerator.h"
#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"

#include "Mesh.h"
#include "SystemShowFunctionSpace.h"
#include "WriterMsh.h"

#include "Exception.h"
#include "Gmsh.h"

using namespace std;

int main(int argc, char** argv){
  // Const
  const char* scalar = "scalar";
  const char* vector = "vector";

  // Init //
  GmshInitialize(argc, argv);

  // Writer //
  WriterMsh writer;

  // Get Domains //
  Mesh msh(argv[1]);
  GroupOfElement domain = msh.getFromPhysical(7);

  // Get FunctionSpace //
  const unsigned int   order = atoi(argv[3]);
  Basis*               basis;
  FunctionSpace*       fSpace;

  if(!strcmp(argv[2], scalar)){
    // If Scalar
    basis =
      BasisGenerator::generate(domain.get(0).getType(),
                               0, order, "hierarchical");

    fSpace = new FunctionSpaceScalar(domain, *basis);
  }

  else if(!strcmp(argv[2], vector)){
    // If Vector
    basis =
      BasisGenerator::generate(domain.get(0).getType(),
                               1, order, "hierarchical");

    fSpace = new FunctionSpaceVector(domain, *basis);
  }

  else
    throw Exception("Unknown FunctionSpace type: %s",
                    argv[2]);

  // Compute all Basis //
  const unsigned int nDof = fSpace->dofNumber();

  for(unsigned int i = 0; i < nDof; i++){
    SystemShowFunctionSpace sys(*fSpace, i);

    sys.assemble();
    sys.solve();

    // View
    stringstream stream;
    stream << "functionSpace_basis" << i;

    writer.setValues(sys);
    writer.write(stream.str());
  }

  // Return //
  GmshFinalize();

  delete fSpace;
  delete basis;

  return 0;
}
