cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose build type")
else(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose build type")
endif(DEFINED CMAKE_BUILD_TYPE)

project(dgplate CXX)

add_subdirectory(../.. "${CMAKE_CURRENT_BINARY_DIR}/gmsh")

include_directories(. ../../Common ../../Numeric ../../Geo ../../Mesh 
   ../../Solver ../../Post ../../Plugin ../../Graphics 
   ${GMSH_EXTERNAL_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/gmsh/Common src ../../contrib/arc ../../contrib/gmm ../dg)

set(SRC
  src/mainDG.cpp
  src/DgC0PlateDof.cpp
  src/IPState.cpp
  src/DgC0PlateSolver.cpp
  src/DgC0PlateTerms.cpp
  src/GModelWithInterface.cpp
  src/LinearElasticShellHookeTensor.cpp
  src/LocalBasis.cpp
  src/MInterfaceElement.cpp
  src/reduction.cpp
  ../dg/dgGroupOfElements.cpp
  ../dg/dgDofContainer.cpp
  ../dg/dgAlgorithm.cpp
  ../dg/dgDataCacheMap.cpp
)

if(HAVE_MPI)
  cmake_force_cxx_compiler(${MPI_COMPILER} "MPI C++ Compiler")
endif(HAVE_MPI)

add_executable(dgplate ${SRC})
target_link_libraries(dgplate shared)

