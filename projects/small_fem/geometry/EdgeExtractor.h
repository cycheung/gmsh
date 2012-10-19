#ifndef _EDGEEXTRACTOR_H_
#define _EDGEEXTRACTOR_H_

#include <map>

#include "Comparators.h"
#include "MElement.h"
#include "MEdge.h"

/**
   @class EdgeExtractor
   @brief A MEdge Extractor

   This class allows @em Edge (MEdge) Extraction from 
   a collection of Elements (MElement).@n

   It got @em only @em class @em methods, so it is @em not requiered
   to instanciate a EdgeExtractor.
*/

class EdgeExtractor{
 public:
   EdgeExtractor(void);
  ~EdgeExtractor(void);

  static std::map<const MEdge*, unsigned int, EdgeComparator>*
    extract(const std::map<const MElement*, 
	                   unsigned int, 
	                   ElementComparator>& element);
  
 private:
  static MEdge* copy(const MEdge& edge);
};


/**
   @fn EdgeExtractor::EdgeExtractor
   Instantiates a new EdgeExtractor
   
   @note
   EdgeExtractor got @em only @em class @em methods, 
   so it is @em not requiered to instanciate it.
   **

   @fn EdgeExtractor::~EdgeExtractor
   Deletes this EdgeExtractor
   **   

   @fn EdgeExtractor::extract
   @param element A map with MElement%s
   @return Returns a map with the MEdge%s in
   the given MElement%s (the mapped values are set to @em zero)
   **
 */

#endif
