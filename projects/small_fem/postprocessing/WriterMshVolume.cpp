#include "WriterMsh.h"

using namespace std;

void WriterMsh::writeVolumeValuesHeader(const string name) const{
  *out << "$ElementData"       << endl
       << "1"                  << endl  // 1 string tag
       << "\"" << name << "\"" << endl  // (name)
       << "1"                  << endl  // 1 real tag
       << "0"                  << endl  // (time value)
       << "3"                  << endl  // 3 integer tag
       << "0"                  << endl  // (time step index)
       << "1"                  << endl  // (number of field -- scalar)
       <<  E                   << endl; // (number of element)
}

void WriterMsh::writeVolumeValues(void) const{
  for(unsigned int i = 0; i < E; i++){
    *out << (*element)[i]->getNum() << " "
         << (*nodalScalarValue)[i]  << " "
         << endl;
  }
}

void WriterMsh::writeVolumeValuesFooter(void) const{
  *out << "$EndElementData" << endl;
}
