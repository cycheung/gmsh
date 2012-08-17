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

   It got @em only @em class @em methods, so it is @em not requiered
   to instanciate a ElementExtractor.
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


/**
   @fn ElementExtractor::ElementExtractor
   Instantiates a new ElementExtractor
   
   @note
   ElementExtractor got @em only @em class @em methods, 
   so it is @em not requiered to instanciate it.
   **

   @fn ElementExtractor::~ElementExtractor
   Deletes this ElementExtractor
   **   

   @fn ElementExtractor::extract
   @param entity A vector of GEntity
   @return Returns an std::pair with:
   @li The first field containing a map with the MElement%s in
   the given entities (the mapped values are set to @em zero)
   @li The second field containing a multimap with the MElement%s
   of the first field, and the @em physicals 
   (see <a href="http://www.geuz.org/gmsh">gmsh</a> documentation) of
   the MElement%s
   **
 */

#endif
