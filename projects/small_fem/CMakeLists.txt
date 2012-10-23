cmake_minimum_required(VERSION 2.6)
project(SmallFem)

## Add Source Macro
macro(add_sources dir files)
  foreach(file ${files})
    list(APPEND list ${dir}/${file})
  endforeach(file)

  set(sources ${sources} ${list} PARENT_SCOPE)
endmacro(add_sources)

## Compatibility with gmsh (TO BE REMOVED !!!)
##############################################
macro(add_sources_in_gmsh dir files)
  foreach(file ${files})
    list(APPEND list ../../${dir}/${file})
  endforeach(file)

  set(sources ${sources} ${list} PARENT_SCOPE)
endmacro(add_sources_in_gmsh)
##############################################

## Look for Packages
include(FindPackageHandleStandardArgs)

## Include gmsh
add_subdirectory(../.. "${CMAKE_CURRENT_BINARY_DIR}/gmsh")

## Add Sources
add_subdirectory(assembler)
add_subdirectory(common)
add_subdirectory(formulation)
add_subdirectory(geometry)
add_subdirectory(postprocessing)
add_subdirectory(solver)
add_subdirectory(simulation)

## TO BE REMOVED !!!
####################
add_subdirectory(../../FunctionSpace 
  "${CMAKE_CURRENT_BINARY_DIR}/FunctionSpace")
####################

## Include Path
include_directories(
  assembler
  common
  formulation
  geometry
  postprocessing
  solver
)

## Include Path for gmsh
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}/gmsh/Common
  ../../Common
  ../../FunctionSpace
  ../../Numeric
  ../../Geo
  ../../Solver
  ${GMSH_EXTERNAL_INCLUDE_DIRS}
)

## Compiler Flags
add_definitions(-pedantic -Wall -g)

## Build SmallFEM Library
add_library(sf STATIC ${sources})
set_target_properties(sf PROPERTIES LINKER_LANGUAGE CXX)

## Build Simulations
add_executable(laplace simulation/Laplace.cpp) 
add_executable(poisson simulation/Poisson.cpp) 

add_executable(projvect    simulation/ProjectionVector.cpp) 
add_executable(projvectsin simulation/ProjectionVectorSin.cpp) 
add_executable(projscalpow simulation/ProjectionScalarPower.cpp) 
add_executable(projscalsin simulation/ProjectionScalarSin.cpp) 
add_executable(projscalone simulation/ProjectionScalarSinOneTime.cpp) 

add_executable(btest simulation/BasisTest.cpp) 
add_executable(conv  simulation/ConvPoisson.cpp) 

## Link (Note: 'lib' is gmsh static library)
target_link_libraries(laplace sf cblas lapack lib ${GMSH_EXTERNAL_LIBRARIES})
target_link_libraries(poisson sf cblas lapack lib ${GMSH_EXTERNAL_LIBRARIES})

target_link_libraries(projvect    sf cblas lapack lib ${GMSH_EXTERNAL_LIBRARIES})
target_link_libraries(projvectsin sf cblas lapack lib ${GMSH_EXTERNAL_LIBRARIES})
target_link_libraries(projscalpow sf cblas lapack lib ${GMSH_EXTERNAL_LIBRARIES})
target_link_libraries(projscalsin sf cblas lapack lib ${GMSH_EXTERNAL_LIBRARIES})
target_link_libraries(projscalone sf cblas lapack lib ${GMSH_EXTERNAL_LIBRARIES})

target_link_libraries(btest sf cblas lapack lib ${GMSH_EXTERNAL_LIBRARIES})
target_link_libraries(conv  sf cblas lapack lib ${GMSH_EXTERNAL_LIBRARIES})
