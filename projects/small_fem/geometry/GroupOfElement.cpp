#include <algorithm>
#include <sstream>
#include <list>

#include "OrientationSort.h"
#include "GroupOfElement.h"

using namespace std;

GroupOfElement::GroupOfElement
(std::multimap<int, const MElement*>::iterator begin,
 std::multimap<int, const MElement*>::iterator end,
 const Mesh& mesh){

  // Get Element //
  list<const MElement*> lst;

  for(; begin != end; begin++)
    lst.push_back(begin->second);

  // Alloc //
  this->mesh = &mesh;
  element = new vector<const MElement*>(lst.begin(), lst.end());
  orientationStat = NULL;
}

GroupOfElement::~GroupOfElement(void){
  delete element;

  if(orientationStat)
    delete orientationStat;
}

void GroupOfElement::
orientAllElements(const Basis& basis){
  // If already oriented, delete old orientation //
  if(orientationStat)
    delete orientationStat;

  // Sort //
  OrientationSort sortPredicate(basis);
  sort(element->begin(), element->end(), sortPredicate);

  // Get Orientation Stats //
  // Get some Data
  const unsigned int nOrient  = basis.getNOrientation();
  const unsigned int nElement = element->size();

  // Init
  orientationStat = new vector<unsigned int>(nOrient);

  for(unsigned int i = 0; i < nOrient; i++)
    (*orientationStat)[i] = 0;

  // Compute
  for(unsigned int i = 0; i < nElement; i++)
    (*orientationStat)[basis.getOrientation(*(*element)[i])]++;

  // The last line is cool isn't it :-) ?
}

string GroupOfElement::toString(void) const{
  stringstream stream;

  stream << "***********************************************"
	 << endl
	 << "* Group Of Element"
	 << endl
	 << "***********************************************"
	 << endl
         << "*                                             *"
	 << endl
	 << "* This group contains the following Elements: *"
	 << endl;

  for(unsigned int i = 0; i < element->size(); i++)
    stream << "*   -- Element #"
	   << mesh->getGlobalId(*(*element)[i])
	   << endl;

  stream << "*                                             *"
	 << endl
         << "***********************************************"
	 << endl;

  return stream.str();
}
