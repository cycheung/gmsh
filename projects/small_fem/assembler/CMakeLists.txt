set(src
  Dof.cpp
  DofManager.cpp
  GroupOfDof.cpp

  AbstractSystem.cpp
  System.cpp
  SystemInstrumented.cpp
  EigenSystem.cpp
)

add_sources(assembler "${src}")
