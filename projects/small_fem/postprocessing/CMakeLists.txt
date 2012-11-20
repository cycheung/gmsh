set(src
  Writer.cpp
  WriterDummy.cpp
  WriterMsh.cpp
  WriterVector.cpp

  Interpolator.cpp
  PlotBasis.cpp

  Integrator.cpp
)


add_sources(postprocessing "${src}")
