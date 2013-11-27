set(src
#  Writer.cpp
#  WriterDummy.cpp
#  WriterVector.cpp

#  WriterMsh.cpp
#  WriterMshElement.cpp
#  WriterMshInterpolation.cpp
#  WriterMshNodal.cpp
#  WriterMshVolume.cpp

  Interpolator.cpp
  FEMSolution.cpp
  ElementSolution.cpp
)


add_sources(postprocessing "${src}")
