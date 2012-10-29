set(src
  Formulation.cpp
  FormulationLaplace.cpp
  FormulationPoisson.cpp
  FormulationProjectionVector.cpp
  FormulationProjectionScalar.cpp

  EigenFormulation.cpp
  FormulationEigenFrequency.cpp
)


add_sources(formulation "${src}")
