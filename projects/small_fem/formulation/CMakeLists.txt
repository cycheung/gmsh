set(src
  FormulationPoisson.cpp
  FormulationProjectionVectorReal.cpp
  FormulationProjectionVectorComplex.cpp
  FormulationProjectionScalarReal.cpp
  FormulationProjectionScalarComplex.cpp
  FormulationSteadyWaveVector.cpp
  FormulationSteadyWaveVectorSlow.cpp
  FormulationSteadyWaveScalarReal.cpp
  FormulationSteadyWaveScalarComplex.cpp

  FormulationNeumann.cpp
  FormulationEMDA.cpp

  FormulationEigenFrequencyVector.cpp
  FormulationEigenFrequencyScalar.cpp

  Term.cpp
  TermFieldField.cpp
  TermGradGrad.cpp
  TermCurlCurl.cpp
  TermProjectionField.cpp
  TermProjectionGrad.cpp

  Quadrature.cpp
)

add_sources(formulation "${src}")
