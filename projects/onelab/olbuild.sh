#!/bin/bash

#The variable 'GMSH_DIR' must be set to the home dir of gmsh on your system
#e.g. export GMSH_DIR=$HOME/SOLVERS/gmsh   

if [ -z "$GMSH_DIR" ]; then
 echo "The variable 'GMSH_DIR' must be set to the home dir of gmsh on your system"
 echo "e.g. export GMSH_DIR=$HOME/SOLVERS/gmsh"
 exit
fi

export COMMON_DIR=$GMSH_DIR/Common
export ONELAB_DIR=$GMSH_DIR/projects/onelab

echo ""
echo "ONELAB: compile onelab interface"

cd SOLVERS
make
cd ..

./olmake.sh
cd ..


