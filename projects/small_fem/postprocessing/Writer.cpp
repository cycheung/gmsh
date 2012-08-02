#include "Writer.h"

Writer::Writer(void){
}

Writer::~Writer(void){
}

void Writer::setValues(std::vector<double>& value){
  hasValue = true;
  isScalar = true;

  nodalScalarValue = &value;
  nodalVectorValue = NULL;
}

void Writer::setValues(std::vector<fullVector<double> >& value){
  hasValue = true;
  isScalar = false;

  nodalScalarValue = NULL;
  nodalVectorValue = &value;
}
