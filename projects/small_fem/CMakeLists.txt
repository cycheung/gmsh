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
)

## Build
add_definitions(-pedantic -Wall -g)
add_executable(laplace ${sources} simulation/Laplace.cpp) 
add_executable(poisson ${sources} simulation/Poisson.cpp) 

add_executable(projvect ${sources} simulation/ProjectionVector.cpp) 
add_executable(projscal ${sources} simulation/ProjectionScalar.cpp) 

add_executable(btest ${sources} simulation/BasisTest.cpp) 
add_executable(conv  ${sources} simulation/ConvPoisson.cpp) 

## Link (Note: 'lib' is gmsh static library)
target_link_libraries(laplace cblas lapack lib)
target_link_libraries(poisson cblas lapack lib)

target_link_libraries(projvect cblas lapack lib)
target_link_libraries(projscal cblas lapack lib)

target_link_libraries(btest cblas lapack lib)
target_link_libraries(conv  cblas lapack lib)
