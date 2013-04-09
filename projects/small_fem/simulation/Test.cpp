#include <iostream>

#include "TriReferenceSpace.h"

#include "Gmsh.h"

using namespace std;


int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  TriReferenceSpace ref;
  cout << ref.toString() << endl;

  GmshFinalize();
  return 0;
}
