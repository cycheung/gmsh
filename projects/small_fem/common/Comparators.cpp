#include "Comparators.h"

bool FaceComparator::operator()(const MFace* a, const MFace* b) const{
  const size_t sizeA = a->getNumVertices();
  const size_t sizeB = b->getNumVertices();

  // Quad Faces are *bigger* than Tri Face
  if(sizeA < sizeB)
    return true;  // 'a' is a Tri and is smaller than 'b' (a quad)

  if(sizeA > sizeB)
    return false; // 'a' is a Quad and is bigger than 'b' (a tri)

  // Here we got both quad or tri
  // --> use Vertex index

  for(size_t i = 0; i < sizeA; i++){
    if(a->getSortedVertex(i) < b->getSortedVertex(i))
      return true;

    if(a->getSortedVertex(i) > b->getSortedVertex(i))
      return false;
  }

  return false;
}
