#include "Exception.h"
#include "Jacobian.h"

using namespace std;

Jacobian::Jacobian(const GroupOfElement& goe,
                   const fullMatrix<double>& point,
		   unsigned int row){

  // Check //
  if((int)row > point.size1() - 1)
    throw Exception("Jacobian: requesting row %d (starting at zero) of a [%d, %d] matrix",
		    row, point.size1(), point.size2());

  if(point.size2() != 3)
    throw Exception("Jacobian: point matrix is of size [N, %d] instead of [N, 3]",
		    point.size1());

  // Get Elements //
  element = &goe.getAll();

  // Get Point //
  this->point = &point;
  this->row   = row;

  // Set maps //
  jac    = NULL;
  invJac = NULL;
}

Jacobian::~Jacobian(void){
  deleteJac();
  deleteInvertJac();
}

void Jacobian::deleteJac(void){
  // Check if NULL //
  if(!jac)
    return;

  // Delete //
  const map<const MElement*, jac_t*, ElementComparator>::iterator
    end = jac->end();

  map<const MElement*, jac_t*, ElementComparator>::iterator
    it = jac->begin();

  for(; it != end; it++){
    delete it->second->first;
    delete it->second;
  }

  delete jac;

  // Set NULL //
  jac = NULL;
}

void Jacobian::deleteInvertJac(void){
  // Check if NULL //
  if(!invJac)
    return;

  // Delete //
  const map<const MElement*, jac_t*, ElementComparator>::iterator
    end = invJac->end();

  map<const MElement*, jac_t*, ElementComparator>::iterator
    it = invJac->begin();

  for(; it != end; it++){
    delete it->second->first;
    delete it->second;
  }

  delete invJac;

  // Set NULL //
  invJac = NULL;
}

void Jacobian::computeJacobians(void){
  // Is already computed ? //
  if(jac)
    return;

  // Init Jac //
  jac = new map<const MElement*, jac_t*, ElementComparator>;

  // Fill Jac //
  jac_t*              tmp;
  fullMatrix<double>* mJac;
  const unsigned int  nElement = element->size();

  for(unsigned int i = 0; i < nElement; i++){
    // Get Jacobian
    tmp  = new jac_t;
    mJac = new fullMatrix<double>(3, 3);

    // GMSH WARNING:
    //    Memory leak with getJacobian() !!

    tmp->second =  const_cast<MElement*>
      ((*element)[i])->getJacobian((*point)(row, 0),
                                   (*point)(row, 1),
                                   (*point)(row, 2),
                                   *mJac);
    tmp->first = mJac;

    // Insert in Jac
    jac->insert(pair<const MElement*, jac_t*>((*element)[i], tmp));
  }
}

void Jacobian::computeInvertFromJac(void){
  // Init InvJac //
  invJac = new map<const MElement*, jac_t*, ElementComparator>;

  // Fill InvJac //
  const map<const MElement*, jac_t*, ElementComparator>::iterator
    end = jac->end();

  map<const MElement*, jac_t*, ElementComparator>::iterator
    it = jac->begin();

  jac_t*              tmp;
  fullMatrix<double>* mIJac;

  for(; it != end; it++){
    tmp = new jac_t;
    mIJac = new fullMatrix<double>(3, 3);

    it->second->first->invert(*mIJac);

    tmp->first  = mIJac;
    tmp->second = it->second->second;

    invJac->insert(pair<const MElement*, jac_t*>(it->first, tmp));
  }
}

void Jacobian::computeInvertFromScratch(void){
  // Init InvJac //
  invJac = new map<const MElement*, jac_t*, ElementComparator>;

  // Fill InvJac //
  jac_t*              tmp;
  fullMatrix<double>* mJac;
  const unsigned int  nElement = element->size();

  for(unsigned int i = 0; i < nElement; i++){
    // Get Jacobian
    tmp  = new jac_t;
    mJac = new fullMatrix<double>(3, 3);

    // GMSH WARNING:
    //    Memory leak with getJacobian() !!

    tmp->second =  const_cast<MElement*>
      ((*element)[i])->getJacobian((*point)(row, 0),
                                   (*point)(row, 1),
                                   (*point)(row, 2),
                                   *mJac);
    mJac->invertInPlace();

    tmp->first = mJac;

    // Insert in InvJac
    invJac->insert(pair<const MElement*, jac_t*>((*element)[i], tmp));
  }
}

const pair<const fullMatrix<double>*, double>&
  Jacobian::getJacobian(const MElement& element) const{

  if(!jac)
    throw Exception("No Jacobian precomputed");

  map<const MElement*, jac_t*, ElementComparator>::iterator it =
    jac->find(&element);

  if(it == jac->end())
    throw
      Exception("Jacobian of Element %d not found",
                element.getNum());

  return *it->second;
}

const pair<const fullMatrix<double>*, double>&
  Jacobian::getInvertJacobian(const MElement& element) const{

  if(!invJac)
    throw Exception("No Invert Jacobian precomputed");

  map<const MElement*, jac_t*, ElementComparator>::iterator it =
    invJac->find(&element);

  if(it == invJac->end())
    throw
      Exception("Invert Jacobian of Element %d not found",
                element.getNum());

  return *it->second;
}
