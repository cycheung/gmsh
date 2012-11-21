#ifndef _GEOEXTRACTOR_H_
#define _GEOEXTRACTOR_H_

#include <vector>
#include <map>

#include "Comparators.h"
#include "GEntity.h"
#include "MElement.h"
#include "MVertex.h"
#include "MEdge.h"
#include "MFace.h"

/**
   @class GeoExtractor
   @brief Extraction of Geometrical Entities

   This class allows the extraction of @em Geomtrical @em Entities:
   @li Element (MElement) extraction from an Entity (GEntity)
   @li Vertex (MVertex) extraction from a collection of Elements (MElement)
   @li Edge (MEdge) extraction from a collection of Elements (MElement)
   @li Face (MFace) extraction from a collection of Elements (MElement)

   It got @em only @em class @em methods, so it is @em not requiered
   to instanciate a GeoExtractor.
*/

class GeoExtractor{
 public:
   GeoExtractor(void);
  ~GeoExtractor(void);
  
  static std::pair<
    std::map<const MElement*, unsigned int, ElementComparator>*,
    std::multimap<int, const MElement*>*
    >
    extractElement(const std::vector<GEntity*>& entity);

  static std::map<const MVertex*, unsigned int, VertexComparator>*
    extractVertex(const std::map<const MElement*, 
		                 unsigned int, 
		                 ElementComparator>& element);

  static std::map<const MEdge*, unsigned int, EdgeComparator>*
    extractEdge(const std::map<const MElement*, 
	                       unsigned int, 
	                       ElementComparator>& element);

  static std::map<const MFace*, unsigned int, FaceComparator>*
    extractFace(const std::map<const MElement*, 
	                       unsigned int, 
	                       ElementComparator>& element);
  
 private:
  static MEdge* copy(const MEdge& edge);
  static MFace* copy(const MFace& face);
};


/**
   @fn GeoExtractor::GeoExtractor
   Instantiates a new GeoExtractor
   
   @note
   GeoExtractor got @em only @em class @em methods, 
   so it is @em not requiered to instanciate it.
   **

   @fn GeoExtractor::~GeoExtractor
   Deletes this GeoExtractor
   **   

   @fn GeoExtractor::extractElement
   @param entity A vector of GEntity
   @return Returns an std::pair with:
   @li The first field containing a map with the MElement%s in
   the given entities (the mapped values are set to @em zero)
   @li The second field containing a multimap with the MElement%s
   of the first field, and the @em physicals 
   (see <a href="http://www.geuz.org/gmsh">gmsh</a> documentation) of
   the MElement%s
   **

   @fn GeoExtractor::extractVertex
   @param element A map with MElement%s
   @return Returns a map with the MVertices in
   the given MElement%s (the mapped values are set to @em zero)
   **

   @fn GeoExtractor::extractEdge
   @param element A map with MElement%s
   @return Returns a map with the MEdge%s in
   the given MElement%s (the mapped values are set to @em zero)
   **

   @fn GeoExtractor::extractFace
   @param element A map with MElement%s
   @return Returns a map with the MFace%s in
   the given MElement%s (the mapped values are set to @em zero)
 */


#endif
