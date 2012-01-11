
This directory contains an example of the coupling of the elasticity solver of gmsh 'elast' with 'gmsh' through the onelab interface.
There are 2 kinds of solver clients for onelab: encapsulated solvers and interfaced solvers.
Encapsulation is to be preferred when the solver client has a decent input data file parser (e.g. getdp),
of which the code is accessible to modification.
For 'elast' has no parser, it is coupled as an interfaced client.

>> cd $GMSH_DIR/projects/onelab/METAMODELS/BEAM
>> gmsh beam.geo

Click  Menu->Tools->Onelab
Enter 'Client Name'=elast 
Browse for the location $GMSH_DIR/projects/onelab/SOLVERS/elast
Click on 'elast'
Click on 'Compute'

It can also be launched non-interactive
>> $GMSH_DIR/projects/onelab/SOLVERS/elast beam 



