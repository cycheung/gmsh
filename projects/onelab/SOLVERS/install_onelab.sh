#!/bin/bash

#Compile onelab_client

cd ../../..
export GMSH_DIR=`pwd`
cd $GMSH_DIR/projects/onelab/SOLVERS
make
echo 'Solver.Executable0 = "_HOME/projects/onelab/SOLVERS/onelab_client";' > zzz
echo 'Solver.Name0 =  "onelab";' >> zzz
sed -e "s|_HOME|${GMSH_DIR}|" zzz > zzz1

#declare onelab_client Solver0 in gmsh
cd
if [ ! -f .gmshrc ]; then 
  echo "creates .gmshrc"
  echo " " > .gmshrc
fi 

cp -f .gmshrc .gmshrc_sav
cat .gmshrc_sav $GMSH_DIR/projects/onelab/SOLVERS/zzz1 > .gmshrc
cd $GMSH_DIR/projects/onelab/SOLVERS
rm -rf zzz zzz1 


