#include "GeoExtractor.h"

GeoExtractor::GeoExtractor(void){
}

GeoExtractor::~GeoExtractor(void){
}

MEdge* GeoExtractor::copy(const MEdge& edge){
  MVertex* begin = edge.getVertex(0);
  MVertex* end   = edge.getVertex(1);

  return new MEdge(begin, end);
}

MFace* GeoExtractor::copy(const MFace& face){
  const int size = face.getNumVertices();
  std::vector<MVertex*> vertex(size);

  for(int i = 0; i < size; i++)
    vertex[i] = face.getVertex(i);

  return new MFace(vertex);
}
