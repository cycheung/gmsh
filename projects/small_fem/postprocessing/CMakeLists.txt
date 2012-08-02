set(src
#  Interpolator.cpp
#  InterpolatorScalar.cpp
#  InterpolatorVector.cpp
#  InterpolatorNode.cpp
#  InterpolatorEdge.cpp
  Writer.cpp
  WriterMsh.cpp

  Solution.cpp
  PlotBasis.cpp
)


add_sources(postprocessing "${src}")
