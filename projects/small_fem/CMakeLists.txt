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
find_package(Doxygen)

## Include gmsh
add_subdirectory(../.. "${CMAKE_CURRENT_BINARY_DIR}/gmsh")

## Add Sources
add_subdirectory(assembler)
add_subdirectory(common)
add_subdirectory(formulation)
add_subdirectory(geometry)
add_subdirectory(postprocessing)
add_subdirectory(solver)

## TO BE REMOVED !!!
####################
add_subdirectory(../../FunctionSpace 
  "${CMAKE_CURRENT_BINARY_DIR}/FunctionSpace")
####################

## Add Main
list(APPEND sources Main.cpp)

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
)

## Build Documentation
if(DOXYGEN_FOUND)
  add_custom_target(doxygen
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile)

endif(DOXYGEN_FOUND)

## Build
add_definitions(-pedantic -Wall -g)
add_executable(sf ${sources}) 

## Link (Note: 'lib' is gmsh static library)
target_link_libraries(sf cblas lapack lib)