#ifndef _ELEMENTEXTRACTOR_H_
#define _ELEMENTEXTRACTOR_H_

#include <vector>
#include <map>

#include "Comparators.h"
#include "GEntity.h"
#include "MElement.h"

/**
   @class ElementExtractor
   @brief A MElement Extractor

   This class allows @em Element (MElement) Extraction from 
   an Entity (GEntity).@n
*/

class ElementExtractor{
 public:
   ElementExtractor(void);
  ~ElementExtractor(void);
  
  static std::pair<
    std::map<const MElement*, unsigned int, ElementComparator>*,
    std::multimap<int, const MElement*>*
    >
    extract(const std::vector<GEntity*>& entity);
};

#endif
