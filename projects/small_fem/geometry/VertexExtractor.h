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

   It got @em only @em class @em methods, so it is @em not requiered
   to instanciate a VertexExtractor.
*/

class VertexExtractor{
 public:
   VertexExtractor(void);
  ~VertexExtractor(void);

  static std::map<const MVertex*, unsigned int, VertexComparator>*
    extract(const std::map<const MElement*, 
	                   unsigned int, 
	                   ElementComparator>& element);
};


/**
   @fn VertexExtractor::VertexExtractor
   Instantiates a new VertexExtractor
   
   @note
   VertexExtractor got @em only @em class @em methods, 
   so it is @em not requiered to instanciate it.
   **

   @fn VertexExtractor::~VertexExtractor
   Deletes this VertexExtractor
   **   

   @fn VertexExtractor::extract
   @param element A map with MElement%s
   @return Returns a map with the MVertices in
   the given MElement%s (the mapped values are set to @em zero)
   **
 */

#endif
