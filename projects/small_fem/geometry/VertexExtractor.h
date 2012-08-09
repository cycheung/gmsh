#ifndef _VERTEXEXTRACTOR_H_
#define _VERTEXEXTRACTOR_H_

#include <map>

#include "Comparators.h"
#include "MElement.h"
#include "MVertex.h"

/**
   @class VertexExtractor
   @brief A MVertex Extractor

   This class allows @em Vertex (MVertex) Extraction from 
   a collection of Elements (MElement).@n
*/

class VertexExtractor{
 public:
   VertexExtractor(void);
  ~VertexExtractor(void);

  static std::map<const MVertex*, unsigned int, MVertexLessThanNum>*
    extract(const std::map<const MElement*, 
	                   unsigned int, 
	                   ElementComparator>& element);
};

#endif
