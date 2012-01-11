
This directory contains examples of the coupling of 'elmer' with 'gmsh' through the onelab interface.

>> cd $GMSH_DIR/projects/onelab/METAMODELS/CRYO

1. Using the generic 'elmer' onelab interface solver

Non interactive:
>> $GMSH_DIR/projects/onelab/SOLVERS/elmer cryo 

Gmsh gui:
>> gmsh cryo.geo
Click  Menu->Tools->Onelab
Enter 'Client Name'=elmer 
Browse for the location $GMSH_DIR/projects/onelab/SOLVERS/elmer
Click on 'elmer'

2. Running the (more sophisticated) *metamodel*, 
that not only computes the solution files 'solution.pos' and 'tempevol.txt' 
but also makes postprocessing treatment.

Compilation
>> ./olmake.sh 

Non interactive:
>> cryosurgery 

Console mode:
>> ../../loader cryosurgery 

