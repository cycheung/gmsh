set(src
  Writer.cpp
  WriterDummy.cpp
  WriterMsh.cpp
  WriterVector.cpp

  Interpolator.cpp
  PlotBasis.cpp
)


add_sources(postprocessing "${src}")
