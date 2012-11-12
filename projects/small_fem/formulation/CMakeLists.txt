set(src
  Formulation.cpp
  FormulationLaplace.cpp
  FormulationPoisson.cpp
  FormulationProjectionVector.cpp
  FormulationProjectionScalar.cpp
  FormulationSteadyWave.cpp

  EigenFormulation.cpp
  FormulationEigenFrequency.cpp
  FormulationVibration.cpp
)


add_sources(formulation "${src}")
