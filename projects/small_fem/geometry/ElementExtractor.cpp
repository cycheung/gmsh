#include "GeoExtractor.h"

using namespace std;

pair<map<const MElement*, size_t, ElementComparator>*,
     multimap<int, const MElement*>*>
GeoExtractor::extractElement(const vector<GEntity*>& entity){
  // Init //
  map<const MElement*, size_t, ElementComparator>*
    element = new map<const MElement*, size_t, ElementComparator>;

  multimap<int, const MElement*>*
    physical = new multimap<int, const MElement*>;

  // Get Elements //
  const size_t nEntity = entity.size();

  for(size_t i = 0; i < nEntity; i++){
    // Get Mesh Elements
    const size_t nElement = entity[i]->getNumMeshElements();
    vector<MElement*> myElement(nElement);

    for(size_t j = 0; j < nElement; j++)
      myElement[j] = entity[i]->getMeshElement(j);

    // Get Physical
    const vector<int> myPhysical = entity[i]->getPhysicalEntities();
    const size_t nPhysical       = myPhysical.size();

    // Insert Element
    for(size_t j = 0; j < nElement; j++){
      pair<map<const MElement*, size_t, ElementComparator>::iterator,
           bool>

        insert = element->insert
        (pair<const MElement*, size_t>(myElement[j], 0));

      // If Insertion is a success,
      // Insert Physical
      if(insert.second)
        for(size_t k = 0; k < nPhysical; k++)
          physical->insert
            (pair<int, const MElement*>(myPhysical[k], myElement[j]));
    }
  }

  // Return //
  return
    pair<map<const MElement*, size_t, ElementComparator>*,
         multimap<int, const MElement*>*>
  (element, physical);
}
