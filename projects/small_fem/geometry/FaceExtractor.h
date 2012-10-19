#ifndef _FACEEXTRACTOR_H_
#define _FACEEXTRACTOR_H_

#include <map>

#include "Comparators.h"
#include "MElement.h"
#include "MFace.h"

/**
   @class FaceExtractor
   @brief A MFace Extractor

   This class allows @em Face (MFace) Extraction from 
   a collection of Elements (MElement).@n

   It got @em only @em class @em methods, so it is @em not requiered
   to instanciate a FaceExtractor.
*/

class FaceExtractor{
 public:
   FaceExtractor(void);
  ~FaceExtractor(void);

  static std::map<const MFace*, unsigned int, FaceComparator>*
    extract(const std::map<const MElement*, 
	                   unsigned int, 
	                   ElementComparator>& element);
  
 private:
  static MFace* copy(const MFace& face);
};


/**
   @fn FaceExtractor::FaceExtractor
   Instantiates a new FaceExtractor
   
   @note
   FaceExtractor got @em only @em class @em methods, 
   so it is @em not requiered to instanciate it.
   **

   @fn FaceExtractor::~FaceExtractor
   Deletes this FaceExtractor
   **   

   @fn FaceExtractor::extract
   @param element A map with MElement%s
   @return Returns a map with the MFace%s in
   the given MElement%s (the mapped values are set to @em zero)
   **
 */

#endif
