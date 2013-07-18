#include <cmath>

#include "BasisGenerator.h"
#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"

#include "Mesh.h"
#include "SystemShowFunctionSpace.h"
#include "WriterMsh.h"
#include "Interpolator.h"

#include "Exception.h"
#include "Gmsh.h"

using namespace std;

int main(int argc, char** argv){
  // Const
  const char* lagrange = "lagrange";
  const char* scalar   = "scalar";
  const char* vector   = "vector";

  // Init //
  GmshInitialize(argc, argv);

  // Writer //
  WriterMsh writer;

  // Get Domains //
  Mesh msh(argv[1]);
  GroupOfElement domain = msh.getFromPhysical(7);

  // Get FunctionSpace //
  const size_t   order = atoi(argv[3]);
  Basis*         basis;
  FunctionSpace* fSpace;

  if(!strcmp(argv[2], lagrange)){
    // If Lagrange
    basis =
      BasisGenerator::generate(domain.get(0).getType(),
                               0, order, "lagrange");

    fSpace = new FunctionSpaceScalar(domain, *basis);
  }

  else if(!strcmp(argv[2], scalar)){
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
  const size_t nDof = fSpace->dofNumber();

  // View names
  char fileName[1024];
  const int nDec = floor(log10(nDof)) + 1;

  for(size_t i = 0; i < nDof; i++){
    SystemShowFunctionSpace sys(*fSpace, i);

    sys.assemble();
    sys.solve();

    // View
    sprintf(fileName, "functionSpace_basis%0*u", nDec, (unsigned int)(i + 1));

    if(argc == 5){
      // Interpolated View //
      Mesh visuMesh(argv[4]);
      GroupOfElement visu = visuMesh.getFromPhysical(7);

      Interpolator interp(sys, visu);
      interp.write(string(fileName), writer);
    }

    else{
      // Adaptive View //
      writer.setValues(sys);
      writer.write(string(fileName));
    }
  }

  // Return //
  GmshFinalize();

  delete fSpace;
  delete basis;

  return 0;
}
