set(src
  Formulation.cpp
  FormulationLaplace.cpp
  FormulationPoisson.cpp
  FormulationProjectionVector.cpp
  FormulationProjectionScalar.cpp
  FormulationSteadyWaveVector.cpp
  FormulationSteadyWaveVectorSlow.cpp
  FormulationSteadyWaveScalar.cpp
  FormulationEigenFrequency.cpp
  FormulationVibration.cpp

  FormulationElementGradGrad.cpp
  FormulationElementOneFormOneForm.cpp
  FormulationElementCurlCurl.cpp

  Term.cpp
  TermFieldField.cpp
  TermGradGrad.cpp
  TermCurlCurl.cpp
  TermProjectionField.cpp
  TermProjectionGrad.cpp

  Quadrature.cpp
)

add_sources(formulation "${src}")
