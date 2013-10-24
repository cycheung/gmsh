#include <cstdlib>
#include <iostream>
#include <sstream>

#include "Mesh.h"
#include "WriterMsh.h"
#include "Gmsh.h"

#include <vector>
#include "Basis.h"
#include "BasisGenerator.h"

using namespace std;

int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  // Get Mesh //
  cout << "## Reading Mesh" << endl << flush;

  Mesh msh(argv[1]);
  size_t       physical = atoi(argv[2]);
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
  WriterMsh writer;
  stringstream name;

  name << "analyze" << physical;

  cout << "## Writing Results" << endl << flush;

  writer.setDomain(domain);
  writer.setValues(orientation);
  writer.write(name.str(), "volume");

  // Done //
  delete basis;
  GmshFinalize();
  cout << "## Done" << endl << flush;
}
