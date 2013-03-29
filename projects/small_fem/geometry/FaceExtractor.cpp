#include "GeoExtractor.h"

using namespace std;

map<const MFace*, unsigned int, FaceComparator>*
GeoExtractor::extractFace(const map<const MElement*,
				    unsigned int,
				    ElementComparator>& element){
  // Init //
  map<const MFace*, unsigned int, FaceComparator>*
    face = new map<const MFace*, unsigned int, FaceComparator>;

  // Get Faces //
  const map<const MElement*, unsigned int, ElementComparator>::const_iterator
    endE = element.end();

  map<const MElement*, unsigned int, ElementComparator>::const_iterator
    itE = element.begin();

  // Iterate on Elements
  for(; itE != endE; itE++){
    // Get Current Element
    MElement* myElement = const_cast<MElement*>(itE->first);

    // Iterate on Faces
    const unsigned int N = myElement->getNumFaces();

    for(unsigned int i = 0; i < N; i++){
      // Take Current Face
      const MFace myFace = myElement->getFace(i);

      // Make a copy (on heap)
      MFace* faceCopy = copy(myFace);

      // Try to Insert
      pair<map<const MFace*, unsigned int, FaceComparator>::iterator,
	   bool> insert =
	face->insert(pair<const MFace* ,int>(faceCopy, 0));

      // If Insertion is not a success,
      // Delete faceCopy
      if(!insert.second)
	delete faceCopy;
    }
  }

  // Return //
  return face;
}
