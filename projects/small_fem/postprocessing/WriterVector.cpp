#include <fstream>
#include <sstream>

#include "Exception.h"
#include "WriterVector.h"

using namespace std;

WriterVector::WriterVector(void){
}

WriterVector::~WriterVector(void){
}

void WriterVector::write(const string name) const{
  // Check if Nodal Value //
  if(!isNodal)
    throw Exception("WriterVector cannot write non Nodal Values");

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

void WriterVector::write(const string name,
                         const string type) const{
  write(name);
}

void WriterVector::write(ostream& stream) const{
  if(isScalar){
    unsigned int size = scalarValue->size();

    for(unsigned int i = 0; i < size; i++)
      stream << i << ": "
	     << scalarValue->at(i)
	     << endl;
  }

  else{
    unsigned int size = vectorValue->size();

    for(unsigned int i = 0; i < size; i++)
      stream << i << ": "
	     << "[ "
	     << (vectorValue->at(i))(0) << " "
	     << (vectorValue->at(i))(1) << " "
	     << (vectorValue->at(i))(2) << " "
	     << " ]"
	     << endl;
  }
}
