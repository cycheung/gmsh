#include <sstream>

#include "Exception.h"
#include "BasisGenerator.h"
#include "WriterMsh.h"

using namespace std;

WriterMsh::WriterMsh(void){
}

WriterMsh::~WriterMsh(void){
}

void WriterMsh::write(const string name) const{
  // Check Domain //
  if(!hasDomain)
    throw Exception("WriterMsh: No Domain has been given !");

  // Open File //
  stringstream fileName;
  fileName << name << ".msh";

  out = new ofstream;
  out->open(fileName.str().c_str());

  // Write Mesh //
  writeHeader();
  writeNodes();
  writeElements();

  // If We have Data to write --> write it :-) //
  if(hasValue){

  // Build InterpolationScheme if needed
    if(!isNodal){
      lBasis = static_cast<BasisLagrange*>
        (BasisGenerator::generate((*element)[0]->getType(),
                                  0,
                                  fs->getBasis(0).getOrder(),
                                  "lagrange")
         );

      writeInterpolationScheme();
    }

    // Header
    writeNodalValuesHeader(name);

    // If NO interpolation needed
    if(isNodal)
      writeNodalValuesFromNode();

    // If interpolation needed
    else
      writeNodalValuesFromSol();

    // Footer
    writeNodalValuesFooter();
  }

  // Close All & Clean //
  out->close();
  delete out;

  if(!isNodal)
    delete lBasis;
}



