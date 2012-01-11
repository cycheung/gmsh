
This directory contains examples of the coupling of 'getdp' with 'gmsh' through the onelab interface.

>> cd $GMSH_DIR/projects/onelab/METAMODELS/CORE

1. Using the generic 'getdp' onelab interface solver

Non interactive:
>> $GMSH_DIR/projects/onelab/SOLVERS/getdp test 

Gmsh gui:
>> gmsh test.geo
Click  Menu->Tools->Onelab
Enter 'Client Name'=getdp 
Browse for the location $GMSH_DIR/projects/onelab/SOLVERS/getdp
Click on 'getdp'

2. Running the *metamodel*, 

Compilation
>> ./olmake.sh 

Non interactive:
>> magcore 

Console mode:
>> ../../loader magcore 




