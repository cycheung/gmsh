# This simple metamodel uses only the post-processor capabilities of gmsh.
# The Physics behind the model is a 2D quasi-static variational 
# ferromagnetic material model. 
# The metamodel successively minimizes a functional depending on 
# the magnetic field (Hx,Hx) and the magnetisation (Jpx, Jpy)
# to determine a new value for Jpx.

# Declaration of the ONELAB parameters of the metamodel:
Hx.number(0, Fields/,"Hx");
Hz.number(0, Fields/,"Hz");
Jpx.number(0, Fields/,"Jpx");
Jpz.number(0, Fields/,"Jpz");
# The definition consists in a Value, a Path and a Label. 
# Take a look to the ONELAB window to see how those elements appear.

# Within a file, the path need not be each time specified, 
# i.e. Hx and Fields/Hx are synonymous in this file. 

# Parameter describing the level of dissipation:
Khi.number(0.12, Parameters/,"Hysteresis coefficient",0:0.5:0.1);
# The 4th argument defines a range (min, max, step) for the parameter
# Range can be seen in the ONELAB window by clicking on the ":" symbol. 

# RadioButton parameters are used for on/off or yes/no options
# They appear as a tickbox in the ONELAB window
DISPLAY.radioButton(1, Options/,"Display functional at each step");


# Gmsh can succcessively call the metamodel in a loop.
# The loop runs over the "list of choices" of the parameter
LOOP.radioButton(1, Options/,"Loop on Hx");
Hx.setAttribute(Loop,OL.get(LOOP));
Hx.resetChoices();
OL.iftrue(DISPLAY)
Hx.addChoices( 0, 0.25, 0.5, 0.75, 1,
    0.76, 0.51, 0.26, 0.01, -0.25, -0.5, -0.75,
    -0.51, -0.01);
OL.else
Hx.addChoices( 0, 0.25, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 
    1.76, 1.51, 1.26, 1.01, 0.76, 0.51, 0.26, 0.01,
    -0.25, -0.5, -0.75, -1, -1.25, -1.50,
    -1.26, -1.01, -0.76, -0.51, -0.26, 0.02);
OL.endif
# There is a shorter list (14 items) when the DISPLAY mode is on
# and a longer list (29 items) when the DISPLAY mode is off.

# One ensures the list of choices of Jpx is reset at the end of the loop
# in case the loop would be run several time. 
OL.if( OL.get(Fields/Hx) == OL.get(Fields/Hx,choices.begin()) )
  Jpx.resetChoices();
OL.endif

# Gmsh can also display on-the-fly graphs of ONELAB variables 
Hx.setAttribute(Graph,10);
Jpx.setAttribute(Graph,01);

# A metamodel is a list of operation:
# pre-processing, computing, post-treatment of the computed data, etc.
# Each operation acts as a client of ONELAB, 
# requesting and sending information to the server.
# In this case, there is only one client, which is the postprocessor of gmsh
# running (non interactively, see the "-" argument) a script called "script"
# This scipt comes with an instrumented version "script.ol"
# (take a look at what happens in this file)
# that is converted by onlob into a valid postprocessing script for gmsh:
Post.register (interfaced, gmsh);
Post.in ( script.ol );
Post.args( script -);
Post.out( minimum.txt );
OL.iftrue(DISPLAY)
  Post.merge(Functional.pos);
OL.else
  Post.merge();
OL.endif
# The postprocessing scripts generates a result file "minimum.txt".
# Relevant info is extracted from the file by ONELAB 
# and transmitted to the server:
Post.up(minimum.txt,1,5,Fields/Jpx, 
        minimum.txt,1,7,Fields/Jpz);


#OL.if( OL.get(Hx) == OL.get(Hx,choices.rbegin()))
#OL.endif