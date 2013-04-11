#include <sstream>

#include "Exception.h"
#include "BasisGenerator.h"
#include "WriterMsh.h"

using namespace std;

WriterMsh::WriterMsh(void){
  lBasis = NULL;
}

WriterMsh::~WriterMsh(void){
}

void WriterMsh::write(const string name) const{
  write(name, "vertex");
}

void WriterMsh::write(const string name,
                      const string type) const{

  if(type.compare("vertex") == 0)
    writeAsNodalValues(name);

  else if(type.compare("volume") == 0)
    writeAsVolumeValues(name);

  else
    throw Exception("WriterMsh: Type %s is not defined!",
                    name.c_str());
}

void WriterMsh::writeInit(const std::string name) const{
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
}

void WriterMsh::writeFinalize(void) const{
  // Close & Clean //
  out->close();
  delete out;

  if(lBasis)
    delete lBasis;
}

void WriterMsh::writeAsNodalValues(const std::string name) const{
  // Init //
  writeInit(name);

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

  // Finalize //
  writeFinalize();
}

void WriterMsh::writeAsVolumeValues(const std::string name) const{
  // Init //
  writeInit(name);

  // If We have Data to write --> write it :-) //
  if(hasValue){
    writeVolumeValuesHeader(name);
    writeVolumeValues();
    writeVolumeValuesFooter();
  }

  // Finalize //
  writeFinalize();
}
