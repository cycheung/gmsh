#include <sstream>
#include "GroupOfJacobian.h"

using namespace std;

GroupOfJacobian::GroupOfJacobian(const GroupOfElement& goe,
                                 const fullMatrix<double>& point,
                                 const string type){
  // Get Elements //
  const unsigned int size = goe.getNumber();

  const vector<const MElement*>& element = goe.getAll();
  this->goe = &goe;

  // Alloc & Populate //
  jac = new Jacobian*[size];

  #pragma omp parallel for
  for(unsigned int i = 0; i < size; i++)
    jac[i] = new Jacobian(*element[i], point, type);
}

GroupOfJacobian::~GroupOfJacobian(void){
  const unsigned int size = goe->getNumber();

  for(unsigned int i = 0; i < size; i++)
    delete jac[i];

  delete[] jac;
}

string GroupOfJacobian::toString(void) const{
  stringstream stream;

  for(unsigned int i = 0; i < goe->getNumber(); i++)
    stream << jac[i]->toString() << endl;

  return stream.str();
}
