#include <cstdio>
#include "DofFixedException.h"

DofFixedException::DofFixedException(const Dof& dof, double value){
  this->dof   = &dof;
  this->value = value;

  string = new char[1024];
  sprintf(string, "Dof %s has a fixed value (%g)",
          dof.toString().c_str(), value);
}

DofFixedException::~DofFixedException(void) throw(){
  delete string;
}

const char* DofFixedException::what(void) const throw(){
  return string;
}

const Dof& DofFixedException::getDof(void) const throw(){
  return *dof;
}

double DofFixedException::getValue(void) const throw(){
  return value;
}
