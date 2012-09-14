Installation of ONELAB

$> export GMSH_DIR=/path/to/your/gmsh/directory
$> cd $GMSH_DIR/projects/onelab/SOLVERS
$> make

The executable $GMSH_DIR/projects/onelab/SOLVERS/onelab_client
should now have been created. It will be used as an external solver by gmsh.


##### A fist example
Let see now how onelab all works with the following tutorial example:

$> cd $GMSH_DIR/projects/onelab/METAMODELS/TUTORIAL
$> gmsh circle.geo -s

In the gmsh menu,
- click on "Tools/Onelab"
- click on the "gear" symbol
- click on "Add new client"
Enter a display name for the external solver (e.g. onelab) and then browse 
the the "onelab_client" executable.

One last setting:
- click on the "gear" symbol again 
- untick "Hide new views"
Now you are ready to interactively exploit the model.
A good first step is to simply click on "Compute". 


### A more involved example with ElmerSolver

If you have ElmerSolver and gnuplot or matlab installed on your system,
you may also take a look to this example:

$> cd $GMSH_DIR/projects/onelab/METAMODELS/LASER
$> gmsh lab_peau.geo -s

