cmake_minimum_required(VERSION 2.6)
project(SmallFem)

## Add Source Macro ##
######################
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

macro(append_gmsh_src dir files)
  foreach(file ${files})
    list(APPEND list ../../${dir}/${file})
  endforeach(file)

  set(sources ${sources} ${list} PARENT_SCOPE)
endmacro(append_gmsh_src)
##############################################

## Look for Packages ##
#######################
include(FindPackageHandleStandardArgs)
include(ExternalProject)

## OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)

## BLAS
if(BLAS_LIBRARIES)
else(BLAS_LIBRARIES)
  find_package(BLAS)
endif(BLAS_LIBRARIES)

## LAPACK
if(LAPACK_LIBRARIES)
else(LAPACK_LIBRARIES)
  find_package(LAPACK)
endif(LAPACK_LIBRARIES)

## Options ##
#############
set(BLAS_LIB   "${BLAS_LIBRARIES}"   CACHE STRING "BLAS libraries path")
set(LAPACK_LIB "${LAPACK_LIBRARIES}" CACHE STRING "LAPACK libraries path")

set(BLAS_LIBRARIES   "${BLAS_LIB}")
set(LAPACK_LIBRARIES "${LAPACK_LIB}")

## Include Gmsh ##
##################
set(BLAS_LAPACK_LIBRARIES "${BLAS_LIBRARIES};${LAPACK_LIBRARIES}")
add_subdirectory(../.. "${CMAKE_CURRENT_BINARY_DIR}/gmsh")

## External Packages ##
#######################
## Linker FLAG for MUMPS
exec_program(${CMAKE_SOURCE_DIR}/script/linker.sh
  ARGS ${BLAS_LIBRARIES} OUTPUT_VARIABLE BLAS_LINKER_FLAGS)

## Install MUMPS
ExternalProject_Add(mumps
  URL http://graal.ens-lyon.fr/MUMPS/MUMPS_4.10.0.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/mumps
  CONFIGURE_COMMAND ""
  BUILD_IN_SOURCE 1
  BUILD_COMMAND make all
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(mumps source_dir)

ExternalProject_Add_Step(mumps setBLAS
  COMMAND ${CMAKE_SOURCE_DIR}/script/replace.sh LIBBLAS "-L${BLAS_LIBRARIES} ${BLAS_LINKER_FLAGS}" ${source_dir}/Make.inc/Makefile.gfortran.SEQ > ${source_dir}/Makefile.inc
  DEPENDEES configure
  DEPENDERS build
)

## MUMPS Include & Libraries
set(MUMPS_INCLUDE ${source_dir}/include)

set(MUMPS_LIB
  ${source_dir}/lib/libcmumps.a
  ${source_dir}/lib/libdmumps.a
  ${source_dir}/lib/libmumps_common.a
  ${source_dir}/lib/libpord.a
  ${source_dir}/lib/libsmumps.a
  ${source_dir}/lib/libzmumps.a
  ${source_dir}/libseq/libmpiseq.a
)

## Add Sources ##
#################
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

## Include Path ##
##################
include_directories(
  ## SmallFEM
  assembler
  common
  formulation
  geometry
  postprocessing
  solver

  ## Gmsh
  ${CMAKE_CURRENT_BINARY_DIR}/gmsh/Common
  ../../Common
  ../../FunctionSpace
  ../../Numeric
  ../../Geo
  ../../Mesh
  ../../Solver

  ## MUMPS
  ${MUMPS_INCLUDE}
)

## Compiler Flags ##
####################
add_definitions(-pedantic -Wall -Wno-long-long) # -Wno-long-long for PETSc 3.3
add_definitions(-O3 -march=native)
add_definitions(-g)

## Build SmallFEM Library ##
############################
## Build
add_library(sf STATIC ${sources})
set_target_properties(sf PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(sf
  ${MUMPS_LIB} lib ${BLAS_LIBRARY} ${LAPACK_LIBRARIES} -lrt)

## Make sure of dependencies
add_dependencies(sf mumps)
add_dependencies(sf lib)

## Build Simulations ##
#######################
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

add_executable(matrix simulation/ElementMatrix.cpp)

## Link ##
##########
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

target_link_libraries(matrix sf)
