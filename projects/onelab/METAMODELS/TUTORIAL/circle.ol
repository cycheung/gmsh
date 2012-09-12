# This simple metamodel uses only gmsh 

# radioButton parameters are used for on/off or yes/no options
DISPLAY.radioButton(0, Options/,"Display functional at each step");

Khi.number(0.12, Parameters/,"Hysteresis coefficient");

# The Fields/Hx parameter is given the value 0

Hx.number(0, Fields/,"Hx"); 

# A list of choices, on the other hand, is used to loop over Fields/Hx.
# There is a shorter list (14 items) when DISPLAY mode is on
# and a longer list (29 items) when the DISPLAY mode is off.
OL.iftrue(DISPLAY)
Hx.addChoices( 0, 0.25, 0.5, 0.75, 1,
               0.76, 0.51, 0.26, 0.01, -0.25, -0.5, -0.75,
	      -0.51, -0.01);
Khi.setValue(0.3,Parameters/);
OL.else
Hx.addChoices( 0, 0.25, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 
               1.76, 1.51, 1.26, 1.01, 0.76, 0.51, 0.26, 0.01,
              -0.25, -0.5, -0.75, -1, -1.25, -1.50,
	      -1.26, -1.01, -0.76, -0.51, -0.26, 0.02);
OL.endif

# other parameters of the metamodel

Hz.number(0, Fields/,"Hz");
Jpx.number(0, Fields/,"Jpx");
Jpz.number(0, Fields/,"Jpz");

# The metamodel has 3 cleints:
# Client 1: Gmsh 

Post.register (interfaced, gmsh);
Post.in ( script.ol );
Post.args( script -);
Post.out( minimum.txt );
OL.iftrue(DISPLAY)
  Post.merge(Functional.pos);
OL.else
  Post.merge();
OL.endif
Post.up(minimum.txt,1,5,Fields/Jpx, 
        minimum.txt,1,7,Fields/Jpz);

# Client 2:

Shell.register(interfaced, sh);
Shell.in( results.sh.ol );
Shell.args( results.sh ); 

OL.if( OL.get(Fields/Hx) == OL.get(Fields/Hx,choices.begin()) )
   Shell.out(results.txt);
OL.else
   Shell.out();
OL.endif

# Remark: Would you define Client 2 like this:
# Shell.register(interfaced, echo);
# Shell.args("OL.get(Fields/Hx) OL.get(Fields/Jpx)" >> results.txt);
# Hx and Jpx would be evaluated at the definition of the client
# and not at run-time. 

# Client 3:

OL.if( OL.get(Fields/Hx) == OL.get(Fields/Hx,choices.rbegin()))
  Gnuplot.register(interfaced, gnuplot);
  Gnuplot.in(results.txt, plot.plt);
  Gnuplot.args(plot.plt);
OL.endif
