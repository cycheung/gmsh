#include <fstream>
#include <sstream>
#include "WriterVector.h"

using namespace std;

WriterVector::WriterVector(void){
}

WriterVector::~WriterVector(void){
}

void WriterVector::write(const string name) const{
  // Check for stdout special case //
  if(name == string("stdout")){
    write(cout);
  }

  // Write to file //
  else{
    stringstream fileName; 
    fileName << name << ".raw";
  
    ofstream out;
    out.open(fileName.str().c_str());
  
    write(out);

    out.close();
  }
}

void WriterVector::write(ostream& stream) const{
  if(isScalar){
    unsigned int size = nodalScalarValue->size();

    for(unsigned int i = 0; i < size; i++)
      stream << i << ": " 
	     << nodalScalarValue->at(i) 
	     << endl;
  }

  else{
    unsigned int size = nodalVectorValue->size();

    for(unsigned int i = 0; i < size; i++)
      stream << i << ": " 
	     << "[ " 
	     << (nodalVectorValue->at(i))(0) << " "
	     << (nodalVectorValue->at(i))(1) << " "
	     << (nodalVectorValue->at(i))(2) << " "
	     << " ]" 
	     << endl;
  }
}
