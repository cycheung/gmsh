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
*/

class EdgeExtractor{
 public:
   EdgeExtractor(void);
  ~EdgeExtractor(void);

  static std::pair<
    std::map<const MEdge*, unsigned int, EdgeComparator>*, 
    std::map<const MEdge*, int, OrientedEdgeComparator>*
    > 
    extract(const std::map<const MElement*, 
	                   unsigned int, 
	                   ElementComparator>& element);

 private:
  static MEdge* copy(const MEdge& edge);
  static MEdge* invert(const MEdge& edge);
};

#endif
