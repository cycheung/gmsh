
This directory contains an example of the coupling of the elasticity solver of gmsh 'elast' with 'gmsh' through the onelab interface.
There are 2 kinds of solver clients for onelab: encapsulated solvers and interfaced solvers.
Encapsulation can be chosen when the solver client has a decent input data file parser,
of which the code is accessible to modification.
In the case of 'elast', there is no parser.
It is then coupled as an interfaced client.

>> cd $GMSH_DIR/projects/onelab/METAMODELS/BEAM

The elast-gmsh interface client is 

$GMSH_DIR/projects/onelab/SOLVERS/elast 

Gmsh gui:
>> gmsh beam.geo
Click  Menu->Tools->Onelab
Enter 'Client Name'=elast 
Browse for the location $GMSH_DIR/projects/onelab/SOLVERS/elast
Click on 'elast'
Click on 'Compute'


It can also be launched non-interactive
>> $GMSH_DIR/projects/onelab/SOLVERS/elast beam 
