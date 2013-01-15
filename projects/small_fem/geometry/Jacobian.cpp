#include "Exception.h"
#include "Jacobian.h"
#include <cstdio>

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
  const vector<const MElement*>&  element = goe.getAll();  
  const unsigned int             nElement = element.size(); 

  // Init Jac //
  jac = new map<const MElement*, jac_t*, ElementComparator>;

  // Fill Jac //
  jac_t* tmp;

  for(unsigned int i = 0; i < nElement; i++){
    // Get Jacobian
    tmp        = new jac_t;
    tmp->first = new fullMatrix<double>(3, 3); 
    
    tmp->second =  const_cast<MElement*>
      (element[i])->getJacobian(point(row, 0), 
				point(row, 1), 
				point(row, 2), 
				*tmp->first);
    // GMSH WARNING:
    //    Memory leak with getJacobian() !!

    // Insert in Jac
    jac->insert(pair<const MElement*, jac_t*>(element[i], tmp));
  }  
}

Jacobian::~Jacobian(void){
 // Delete Jac //
  const map<const MElement*, jac_t*, ElementComparator>::iterator 
    endJ = jac->end();

  map<const MElement*, jac_t*, ElementComparator>::iterator
    itJ = jac->begin();
 
  for(; itJ != endJ; itJ++){
    delete itJ->second->first;
    delete itJ->second;    
  }

  delete jac;
}
