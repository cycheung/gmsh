#include <iostream>

#include "SmallFem.h"

#include "Timer.h"
#include "TriLagrangeBasis.h"
#include "LineReferenceSpace.h"
#include "TriReferenceSpace.h"
#include "QuadReferenceSpace.h"
#include "TetReferenceSpace.h"
#include "HexReferenceSpace.h"
#include "PyrReferenceSpace.h"
#include "PriReferenceSpace.h"

#include "LineNodeBasis.h"
#include "LineEdgeBasis.h"
#include "LineNedelecBasis.h"
#include "TriNodeBasis.h"
#include "QuadNedelecBasis.h"

#include "Mesh.h"
#include "fullMatrix.h"
#include "GroupOfJacobian.h"

#include "PermutationTree.h"

#include "SolverMatrix.h"
#include "SolverVector.h"
#include "SolverMUMPS.h"

using namespace std;

int main(int argc, char** argv){
  SmallFem::Keywords("-msh,-o");
  SmallFem::Initialize(argc, argv);

  Options& option = SmallFem::getOptions();
  cout << option.toString() << endl;

  SmallFem::Finalize();
}
