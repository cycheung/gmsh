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

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

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
  ../../Mesh
  ../../Solver
  ${GMSH_EXTERNAL_INCLUDE_DIRS}
)

## Include Path for MUMPS
include_directories(
  contrib
)

## Library Path for MUMPS
set(mumps
  ${CMAKE_SOURCE_DIR}/contrib/libcmumps.a
  ${CMAKE_SOURCE_DIR}/contrib/libdmumps.a
  ${CMAKE_SOURCE_DIR}/contrib/libmumps_common.a
  ${CMAKE_SOURCE_DIR}/contrib/libpord.a
  ${CMAKE_SOURCE_DIR}/contrib/libsmumps.a
  ${CMAKE_SOURCE_DIR}/contrib/libzmumps.a
  ${CMAKE_SOURCE_DIR}/contrib/libmpiseq.a
)

## Compiler Flags
add_definitions(-pedantic -Wall -Wno-long-long) # -Wno-long-long for PETSc 3.3
add_definitions(-O3 -march=native)
add_definitions(-g)

## Build SmallFEM Library
add_library(sf STATIC ${sources})
set_target_properties(sf PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(sf ${mumps} lib ${GMSH_EXTERNAL_LIBRARIES} -lrt)

## Build Simulations
add_executable(laplace simulation/Laplace.cpp)
add_executable(poisson simulation/Poisson.cpp)

add_executable(projvectsin simulation/ProjectionVectorSin.cpp)
add_executable(projvectone simulation/ProjectionVectorSinOneTime.cpp)

add_executable(projscalsin simulation/ProjectionScalarSin.cpp)
add_executable(projscalone simulation/ProjectionScalarSinOneTime.cpp)

add_executable(swavev    simulation/SteadyWaveVector.cpp)
add_executable(swaves    simulation/SteadyWaveScalar.cpp)
add_executable(cavity    simulation/Cavity.cpp)
add_executable(vibration simulation/Vibration.cpp)

add_executable(mesh simulation/Mesh.cpp)
add_executable(fs   simulation/ShowFunctionSpace.cpp)
add_executable(test simulation/Test.cpp)

## Link (Note: 'lib' is gmsh static library)
target_link_libraries(laplace sf)
target_link_libraries(poisson sf)

target_link_libraries(projvectsin sf)
target_link_libraries(projvectone sf)

target_link_libraries(projscalsin sf)
target_link_libraries(projscalone sf)

target_link_libraries(swavev    sf)
target_link_libraries(swaves    sf)
target_link_libraries(cavity    sf)
target_link_libraries(vibration sf)

target_link_libraries(mesh sf)
target_link_libraries(fs   sf)
target_link_libraries(test sf)
