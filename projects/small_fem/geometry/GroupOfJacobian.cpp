#include <sstream>
#include "GroupOfJacobian.h"

using namespace std;

GroupOfJacobian::GroupOfJacobian(const GroupOfElement& goe,
                                 const Basis& basis,
                                 const fullMatrix<double>& point,
                                 const string type){
  // Get Elements //
  const size_t size = goe.getNumber();

  const vector<const MElement*>& element = goe.getAll();
  this->goe = &goe;

  // Alloc & Populate //
  jac = new Jacobian*[size];

  // Compute first jacobian outside OpenMP
  // Need to generate GMSH nodalBasis outside a thread
  // so that Vandermond matrix is not inverted
  // in nested threads (limitation of OpenBLAS)
  jac[0] = new Jacobian(*element[0], basis, point, type);

  #pragma omp parallel for
  for(size_t i = 1; i < size; i++)
    jac[i] = new Jacobian(*element[i], basis, point, type);
}

GroupOfJacobian::~GroupOfJacobian(void){
  const size_t size = goe->getNumber();

  for(size_t i = 0; i < size; i++)
    delete jac[i];

  delete[] jac;
}

string GroupOfJacobian::toString(void) const{
  stringstream stream;

  for(size_t i = 0; i < goe->getNumber(); i++)
    stream << jac[i]->toString() << endl;

  return stream.str();
}
