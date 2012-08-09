#include <sstream>
#include <list>

#include "GroupOfElement.h"

using namespace std;

unsigned int GroupOfElement::nextId = 0;

GroupOfElement::GroupOfElement
(std::multimap<int, const MElement*>::iterator begin, 
 std::multimap<int, const MElement*>::iterator end){
  
  // Set ID //
  id = nextId;
  nextId++;
  
  // Get Element //
  list<const MElement*> lst;

  for(; begin != end; begin++)
    lst.push_back(begin->second);

  // Alloc //
  nElement = lst.size();
  element  = 
    new vector<const MElement*>(lst.begin(), lst.end());
}

GroupOfElement::~GroupOfElement(void){
  delete element;
}

string GroupOfElement::toString(void) const{
  stringstream stream;
  
  stream << "*********************************************"    
	 << endl
	 << "* Group Of Element #" << id    
	 << endl
	 << "*********************************************" 
	 << endl << "*" 
	 << endl
	 << "* This group contains "
	 << nElement
	 << " elements" 
	 << endl
	 << "*********************************************" 
	 << endl;
  
  return stream.str();
}
