cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose build type")
else(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose build type")
endif(DEFINED CMAKE_BUILD_TYPE)

# variable for other projects (not include NonLinearMechSolver)
set(CM3APPS TRUE)

if(CMAKE_BUILD_TYPE MATCHES "debug")
  set(INDEBUG TRUE)
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
  set(INDEBUG TRUE)
else()
  set(INDEBUG FALSE)
endif(CMAKE_BUILD_TYPE MATCHES "debug")

if(INDEBUG)
  set(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} -D_DEBUG")
  set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -D_DEBUG")
  set(LIBDGSHELLNAME "dgshellpyDebug")
else()
#  set(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} -g -pg")
#  set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -g -pg")
  set(LIBDGSHELLNAME "dgshellpy")
endif(INDEBUG)

project(cm3apps)
# includes or not the other projects (default no project)
option(ENABLE_DGSHELL "dgshell" OFF)
option(ENABLE_DG3D "dG3D" OFF)
option(ENABLE_MSCH "msch" OFF)

#  NonLinearSolver is mandatory
add_subdirectory(../NonLinearSolver "${CMAKE_CURRENT_BINARY_DIR}/NonLinearSolver")
include_directories(. ../../Common ../../Numeric ../../Geo ../../Mesh
   ../../Solver ../../Post ../../Plugin ../../Graphics ../NonLinearSolver/internalPoints ../NonLinearSolver/nlsolver ../NonLinearSolver ../NonLinearSolver/Interface ../NonLinearSolver/contact ../NonLinearSolver/materialLaw ../NonLinearSolver/Domain ../NonLinearSolver/BoundaryConditions ../NonLinearSolver/nlTerms ../NonLinearSolver/field ../NonLinearSolver/space ${GMSH_EXTERNAL_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/NonLinearSolver/gmsh/Common src ../../contrib/gmm ${CMAKE_CURRENT_BINARY_DIR}/NonLinearSolver ../NonLinearSolver/periodicBC $ENV{NLSMPIINC} )

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/src *.h)

if(HAVE_MPI)
  cmake_force_cxx_compiler(${MPI_COMPILER} "MPI C++ Compiler")
endif(HAVE_MPI)

if(ENABLE_DGSHELL)
  add_subdirectory(../dgshell "${CMAKE_CURRENT_BINARY_DIR}/dgshell")
endif(ENABLE_DGSHELL)

if(ENABLE_DG3D)
  add_subdirectory(../dG3D "${CMAKE_CURRENT_BINARY_DIR}/dG3D")
endif(ENABLE_DG3D)

if(ENABLE_MSCH)
  add_subdirectory(../msch "${CMAKE_CURRENT_BINARY_DIR}/msch")
endif(ENABLE_MSCH)
