set(src
  Writer.cpp
  WriterDummy.cpp
  WriterVector.cpp

  WriterMsh.cpp
  WriterMshElement.cpp
  WriterMshInterpolation.cpp
  WriterMshNodal.cpp
  WriterMshVolume.cpp

  Interpolator.cpp
  PlotBasis.cpp
)


add_sources(postprocessing "${src}")
