set(src
  Dof.cpp
  DofManager.cpp
  DofFixedException.cpp
  GroupOfDof.cpp

  SystemAbstract.cpp
  System.cpp
  SystemEigen.cpp
  SystemInstrumented.cpp
  SystemShowFunctionSpace.cpp
)

add_sources(assembler "${src}")
