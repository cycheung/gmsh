#!/bin/bash

#The variable 'GMSH_DIR' must be set to the home dir of gmsh on your system
#e.g. export GMSH_DIR=$HOME/SOLVERS/gmsh   

if [ -z "$GMSH_DIR" ]; then
 echo "The variable 'GMSH_DIR' must be set to the home dir of gmsh on your system"
 echo "e.g. export GMSH_DIR=$HOME/SOLVERS/gmsh"
 exit
fi

if [ -z  "$`which gmsh`" ]; then
 echo "The command 'gmsh' is not defined systemwide"
fi
if [ -z  "`which getdp`" ]; then
 echo "The command 'getdp' is not defined systemwide"
fi
if [ -z  "`which elmer`" ]; then
 echo "The command 'elmer' is not defined systemwide"
fi

export COMMON_DIR=$GMSH_DIR/Common
export ONELAB_DIR=$GMSH_DIR/projects/onelab

echo ""
echo "ONELAB: compile the console (non gui) onelab loader"
g++ -I $ONELAB_DIR -I $COMMON_DIR $COMMON_DIR/OS.cpp $ONELAB_DIR/OnelabClients.cpp $ONELAB_DIR/OnelabMessage.cpp $ONELAB_DIR/loader.cpp -o loader

exit 

echo "ONELAB: compile the metamodels"

cd METAMODELS

cd CORE
olmake.sh
cd ..

cd CRYO
olmake.sh
cd ..

cd ..

echo "ONELAB: compile a simple remote client for testing socket communication" 
g++ -I $ONELAB_DIR -I $COMMON_DIR remote.cpp -o remote
