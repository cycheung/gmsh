#include <algorithm>
#include <sstream>
#include <list>

#include "Exception.h"
#include "GroupOfElement.h"

using namespace std;

GroupOfElement::OrientationSort::OrientationSort(const Basis& basis){
  this->basis = &basis;
}

GroupOfElement::OrientationSort::~OrientationSort(void){
}

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
  element.assign(lst.begin(), lst.end());

  orientationStat.clear();
}

GroupOfElement::~GroupOfElement(void){
}

size_t GroupOfElement::getNumber(void) const{
  return element.size();
}

const std::vector<const MElement*>& GroupOfElement::getAll(void) const{
  return element;
}


const Mesh& GroupOfElement::getMesh(void) const{
  return *mesh;
}

void GroupOfElement::
orientAllElements(const Basis& basis){
  // If already oriented, delete old orientation //
  if(orientationStat.size() != 0)
    throw Exception
      ("GroupOfElement already oriented");

  // Sort //
  OrientationSort predicate(basis);
  sort(element.begin(), element.end(), predicate);

  // Get Orientation Stats //
  // Get some Data
  const size_t nOrient  = basis.getReferenceSpace().getNReferenceSpace();
  const size_t nElement = element.size();

  // Init
  orientationStat.resize(nOrient);

  for(size_t i = 0; i < nOrient; i++)
    orientationStat[i] = 0;

  // Compute
  for(size_t i = 0; i < nElement; i++)
    orientationStat[basis.getReferenceSpace().getReferenceSpace(*element[i])]++;

  // The last line is cool isn't it :-) ?
}

const std::vector<size_t>& GroupOfElement::getOrientationStats(void) const{
  if(orientationStat.size() != 0)
    return orientationStat;

  else
    throw Exception("Orientation of GroupOfElement not computed");
}


void GroupOfElement::unoriented(void){
  orientationStat.clear();
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

  for(size_t i = 0; i < element.size(); i++)
    stream << "*   -- Element #"
           << mesh->getGlobalId(*element[i])
           << endl;

  stream << "*                                             *"
         << endl
         << "***********************************************"
         << endl;

  if(orientationStat.size() == 0)
    return stream.str();


  stream << "*                                             *"
         << endl
         << "* This group has the following Orientations:  *"
         << endl;

  for(size_t i = 0; i < orientationStat.size(); i++)
    stream << "*   -- Elements with Orientation " << i << " - "
           << orientationStat[i] << endl;

  stream << "*                                             *"
         << endl
         << "***********************************************"
         << endl;

  return stream.str();
}
