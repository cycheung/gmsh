#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include "Mesh.h"
#include "Basis.h"
#include "BasisGenerator.h"

#include "ElementSolution.h"
#include "SmallFem.h"

using namespace std;

void compute(const Options& option){
  // Get Mesh //
  cout << "## Reading Mesh" << endl << flush;

  Mesh msh(option.getValue("-msh")[1]);
  size_t       physical = atoi(option.getValue("-phys")[1].c_str());
  GroupOfElement domain = msh.getFromPhysical(physical);

  // Select Basis //
  cout << "## Generate Basis" << endl << flush;
  Basis* basis =
    BasisGenerator::generate(domain.get(0).getType(),
                             0, 1, "hierarchical");
  // Get Reference Space //
  const ReferenceSpace& refSpace = basis->getReferenceSpace();

  // Alloc Orientations //
  const size_t size = domain.getNumber();
  vector<double> orientation(size);

  // Analyze //
  cout << "## Analyzing" << endl << flush;

  for(size_t i = 0; i < size; i++)
    orientation[i] = (double)(refSpace.getReferenceSpace(domain.get(i)));

  // Full Mesh //
  cout << "## Full Mesh data" << endl << flush
       << msh.toString()      << endl << flush;

  // Mesh //
  cout << "## Mesh" << endl << flush;
  for(size_t i = 0; i < domain.getNumber(); i++){
    const MElement& e = domain.get(i);
    const size_t    N = e.getNumPrimaryVertices();

    vector<int> v(N);
    for(size_t j = 0; j < N; j++)
      v[j] = e.getVertex(j)->getNum();

    cout << "  -- " << "Element " << e.getNum()
         << ": "    << "[";

    for(size_t j = 0; j < N - 1; j++)
      cout << v[j] << ", ";

    cout << v[N - 1] << "]"
         << ": #" << orientation[i]
         << endl;
  }

  // Printing //
  stringstream name;
  name << "analyze" << physical;

  cout << "## Writing Results" << endl << flush;

  ElementSolution sol;
  sol.addValues(0, 0, domain, orientation);
  sol.write(name.str());

  // Done //
  delete basis;
  cout << "## Done" << endl << flush;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-phys");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
