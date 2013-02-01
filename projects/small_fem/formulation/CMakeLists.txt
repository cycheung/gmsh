set(src
  Formulation.cpp
  FormulationLaplace.cpp
  FormulationPoisson.cpp
  FormulationProjectionVector.cpp
  FormulationProjectionScalar.cpp
  FormulationSteadyWaveVector.cpp
  FormulationSteadyWaveScalar.cpp

  SlowFormulationSteadyWaveVector.cpp

  EigenFormulation.cpp
  FormulationEigenFrequency.cpp
  FormulationVibration.cpp

  Term.cpp
  TermFieldField.cpp
  TermGradGrad.cpp
  TermCurlCurl.cpp
  TermProjectionField.cpp
  TermProjectionGrad.cpp
)


add_sources(formulation "${src}")
