set(src
  FormulationPoisson.cpp
  FormulationProjectionVector.cpp
  FormulationProjectionScalar.cpp
  FormulationSteadyWaveVector.cpp
  FormulationSteadyWaveVectorSlow.cpp
  FormulationSteadyWaveScalar.cpp

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
