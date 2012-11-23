set(src
  Formulation.cpp
  FormulationLaplace.cpp
  FormulationPoisson.cpp
  FormulationProjectionVector.cpp
  FormulationProjectionScalar.cpp
  FormulationSteadyWaveVector.cpp
  FormulationSteadyWaveScalar.cpp

  EigenFormulation.cpp
  FormulationEigenFrequency.cpp
  FormulationVibration.cpp
)


add_sources(formulation "${src}")
