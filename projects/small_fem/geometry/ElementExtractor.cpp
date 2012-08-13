#include "ElementExtractor.h"

using namespace std;

ElementExtractor::ElementExtractor(void){
}


ElementExtractor::~ElementExtractor(void){
}

pair<map<const MElement*, unsigned int, ElementComparator>*,
     multimap<int, const MElement*>*>

ElementExtractor::extract(const vector<GEntity*>& entity){
  // Init //
  map<const MElement*, unsigned int, ElementComparator>* 
    element = new map<const MElement*, unsigned int, ElementComparator>;

  multimap<int, const MElement*>* 
    physical = new multimap<int, const MElement*>;

  // Get Elements //
  const unsigned int nEntity = entity.size();

  for(unsigned int i = 0; i < nEntity; i++){
    // Get Mesh Elements
    const unsigned int nElement = entity[i]->getNumMeshElements();
    vector<MElement*> myElement(nElement);

    for(unsigned int j = 0; j < nElement; j++)
      myElement[j] = entity[i]->getMeshElement(j);

    // Get Physical
    const vector<int> myPhysical = entity[i]->getPhysicalEntities();
    const unsigned int nPhysical = myPhysical.size();

    // Insert Element     
    for(unsigned int j = 0; j < nElement; j++){
      pair<map<const MElement*, unsigned int, ElementComparator>::iterator,
	   bool>
	
	insert = element->insert
	(pair<const MElement*, unsigned int>(myElement[j], 0));
     
      // If Insertion is a success,
      // Insert Physical
      if(insert.second)
	for(unsigned int k = 0; k < nPhysical; k++)
	  physical->insert
	    (pair<int, const MElement*>(myPhysical[k], myElement[j]));  
    }    
  }

  // Return //
  return 
    pair<map<const MElement*, unsigned int, ElementComparator>*,
	 multimap<int, const MElement*>*>
  (element, physical);
}
