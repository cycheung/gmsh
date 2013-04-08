#include <iostream>

#include "TriEdgeBasis.h"

#include "Gmsh.h"

using namespace std;


int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  TriEdgeBasis(atoi(argv[1]));

  GmshFinalize();
  return 0;
}
