#include "Exception.h"
#include "Jacobian.h"

using namespace std;

Jacobian::Jacobian(const GroupOfElement& goe,
                   const fullMatrix<double>& point){

  // Check //
  if(point.size2() != 3)
    throw Exception("Jacobian: point matrix is of size [N, %d] instead of [N, 3]",
		    point.size1());

  // Get Elements //
  element  = &goe.getAll();
  nElement = element->size();

  // Get Point //
  this->point = &point;
  nPoint      = point.size1();

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
    for(unsigned int i = 0; i < nPoint; i++){
      delete (*it->second)[i]->first;
      delete (*it->second)[i];
    }

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
    for(unsigned int i = 0; i < nPoint; i++){
      delete (*it->second)[i]->first;
      delete (*it->second)[i];
    }

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
  jac_pair*           tmpJac_pair;
  jac_t*              tmpJac_t;
  fullMatrix<double>* mJac;

  // Loop on Element
  for(unsigned int i = 0; i < nElement; i++){
    tmpJac_t  = new jac_t(nPoint);

    // Loop on Points
    for(unsigned int j = 0; j < nPoint; j++){
      tmpJac_pair = new jac_pair;
      mJac        = new fullMatrix<double>(3, 3);

      tmpJac_pair->second = const_cast<MElement*>
        ((*element)[i])->getJacobian((*point)(j, 0),
                                     (*point)(j, 1),
                                     (*point)(j, 2),
                                     *mJac);
      tmpJac_pair->first = mJac;
      (*tmpJac_t)[j]     = tmpJac_pair;
    }

    // Insert in Jac
    jac->insert(pair<const MElement*, jac_t*>((*element)[i], tmpJac_t));
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

  jac_pair*           tmpJac_pair;
  jac_t*              tmpJac_t;
  fullMatrix<double>* mIJac;

  // Loop on Elements
  for(; it != end; it++){
    tmpJac_t = new jac_t(nPoint);

    // Loop on Points
    for(unsigned int j = 0; j < nPoint; j++){
      tmpJac_pair = new jac_pair;
      mIJac       = new fullMatrix<double>(3, 3);

      (*it->second)[j]->first->invert(*mIJac);

      tmpJac_pair->first  = mIJac;
      tmpJac_pair->second = (*it->second)[j]->second;

      (*tmpJac_t)[j] = tmpJac_pair;
    }

    // Insert in invJac
    invJac->insert(pair<const MElement*, jac_t*>(it->first, tmpJac_t));
  }
}

void Jacobian::computeInvertFromScratch(void){
  // Init InvJac //
  invJac = new map<const MElement*, jac_t*, ElementComparator>;

  // Fill InvJac //
  jac_pair*           tmpJac_pair;
  jac_t*              tmpJac_t;
  fullMatrix<double>* mIJac;

  // Loop on Element
  for(unsigned int i = 0; i < nElement; i++){
    tmpJac_t  = new jac_t(nPoint);

    // Loop on Points
    for(unsigned int j = 0; j < nPoint; j++){
      tmpJac_pair = new jac_pair;
      mIJac       = new fullMatrix<double>(3, 3);

      tmpJac_pair->second =  const_cast<MElement*>
        ((*element)[i])->getJacobian((*point)(j, 0),
                                     (*point)(j, 1),
                                     (*point)(j, 2),
                                     *mIJac);
      mIJac->invertInPlace();

      tmpJac_pair->first = mIJac;
      (*tmpJac_t)[j]     = tmpJac_pair;
    }

    // Insert in InvJac
    invJac->insert(pair<const MElement*, jac_t*>((*element)[i], tmpJac_t));
  }
}

const vector<const pair<const fullMatrix<double>*, double>*>&
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

const vector<const pair<const fullMatrix<double>*, double>*>&
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
