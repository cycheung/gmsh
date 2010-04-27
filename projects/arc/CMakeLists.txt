cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
option(ENABLE_NETCDF "Enable NETCDF" ON)

# if CMAKE_BUILD_TYPE is specified use it; otherwise set the default
# build type to "RelWithDebInfo" ("-O2 -g" with gcc) prior to calling
# project()
if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose build type")
else(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose build type")
endif(DEFINED CMAKE_BUILD_TYPE)

project(arc CXX)

add_subdirectory(../.. "${CMAKE_CURRENT_BINARY_DIR}/gmsh")


FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)


#list(APPEND EXTERNAL_LIBRARIES -lgsl)

include_directories(. ../../Common ../../Numeric ../../Geo ../../Mesh
   ../../Solver ../../Post ../../Plugin ../../Graphics ../../contrib/gmm ../../contrib/DiscreteIntegration
   ${GMSH_EXTERNAL_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/gmsh/Common ${EXTERNAL_INCLUDES})

#set(SRC
#  mainElasticity.cpp
#)

message(ginclude ${GMSH_EXTERNAL_INCLUDES_DIRS})
set (LINK_LIBRARIES shared ${EXTERNAL_LIBRARIES})

add_executable(mainElasticity EXCLUDE_FROM_ALL mainElasticity.cpp ${CMAKE_CURRENT_BINARY_DIR}/gmsh)
target_link_libraries(mainElasticity ${LINK_LIBRARIES})

add_executable(XFEMInclusion EXCLUDE_FROM_ALL XFEMInclusion.cpp ${CMAKE_CURRENT_BINARY_DIR}/gmsh)
target_link_libraries(XFEMInclusion ${LINK_LIBRARIES})

add_executable(ImportLSImage EXCLUDE_FROM_ALL ImportLSImage.cpp ${CMAKE_CURRENT_BINARY_DIR}/gmsh)
target_link_libraries(ImportLSImage ${LINK_LIBRARIES} ITKIO)

add_executable(ImportLS2dImage EXCLUDE_FROM_ALL ImportLS2dImage.cpp ${CMAKE_CURRENT_BINARY_DIR}/gmsh)
target_link_libraries(ImportLS2dImage ${LINK_LIBRARIES} ITKIO)

add_executable(mainImageSolver EXCLUDE_FROM_ALL mainImageSolver.cpp ${CMAKE_CURRENT_BINARY_DIR}/gmsh)
target_link_libraries(mainImageSolver ${LINK_LIBRARIES} ITKIO)

add_executable(mainBarElement EXCLUDE_FROM_ALL mainBarElement.cpp ${CMAKE_CURRENT_BINARY_DIR}/gmsh)
target_link_libraries(mainBarElement ${LINK_LIBRARIES})

add_executable(mainGui EXCLUDE_FROM_ALL mainGui.cpp ${CMAKE_CURRENT_BINARY_DIR}/gmsh)
target_link_libraries(mainGui ${LINK_LIBRARIES} ITKIO ${FLTK_LIBRARIES})


IF ( FLTK_FOUND )
       FLTK_WRAP_UI(arc arc.fl )
ENDIF ( FLTK_FOUND )




