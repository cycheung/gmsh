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
if [ -z  "`which ElmerGrid`" ]; then
 echo "The command 'ElmerGrid' is not defined systemwide"
fi
if [ -z  "`which ElmerSolver`" ]; then
 echo "The command 'ElmerSolver' is not defined systemwide"
fi

export COMMON_DIR=$GMSH_DIR/Common
export ONELAB_DIR=$GMSH_DIR/projects/onelab

echo ""
echo "ONELAB: compile onelag generic solvers"

cd SOLVERS
./olmake.sh
cd ..

cd METAMODELS
echo ""
echo "ONELAB: compile the metamodels"

cd CORE
./olmake.sh
cd ..

cd CRYO
./olmake.sh
cd ..

cd ..
./olmake.sh

