
# Declaration of the ONELAB parameters of the metamodel:
Hx.number(0, Fields/,"Hx");
Hz.number(0, Fields/,"Hz");
Jpx.number(0, Fields/,"Jpx");
Jpz.number(0, Fields/,"Jpz");
Khi.number(0.12, Parameters/,"Hysteresis coefficient",0:0.5:0.1);
DISPLAY.radioButton(1, Options/,"Display functional at each step");


Mesh.register(encapsulated, gmsh);
Mesh.in( circle.geo);
Mesh.out(circle.msh);
Mesh.args(circle.geo);

# In this metamodel, the only client is the postprocessor of gmsh
# running (non interactively, see the "-" argument) a script called "script"
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

# Within a file, the path need not be each time specified, 
# i.e. Hx and Fields/Hx are synonymous in this file. 

# One ensures the list of choices of Jpx is reset at the beginning of the loop
# in case the loop would be run several time. 
OL.if( OL.get(Fields/Hx) == OL.get(Fields/Hx,choices.begin()) )
  Jpx.resetChoices();
OL.endif

# Gmsh can also display on-the-fly graphs of ONELAB variables 
Hx.setAttribute(Graph,1000);
Jpx.setAttribute(Graph,0100);
