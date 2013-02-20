set(src
  Mapper.cpp
  Mesh.cpp
  GroupOfElement.cpp

  GeoExtractor.cpp
  ElementExtractor.cpp
  VertexExtractor.cpp
  EdgeExtractor.cpp
  FaceExtractor.cpp

  OrientationSort.cpp

  GroupOfJacobian.cpp
  Jacobian.cpp
)

add_sources(geometry "${src}")
