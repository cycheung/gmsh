
# radioButton parameters are used for on/off or yes/no options
DISPLAY.radioButton(0,Options/,"Display functional at each step");

Khi.number(0.12, Parameters/,"Hysteresis coefficient");

# The Fields/Hs parameter is given the default value 0.5
# The list of choices, o, the other hand, is used when looping on Hx.
Hx.number(0, Fields/,"Hx"); 

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

Hz.number(0, Fields/,"Hz");
Jpx.number(0, Fields/,"Jpx");
Jpz.number(0, Fields/,"Jpz");

# Client 1: 

Post.register (interfaced, gmsh);
Post.in ( script.ol );
Post.args( script -);
Post.out( minimum.txt );
OL.iftrue(DISPLAY)
  Post.merge(Functional.pos);
OL.else
  Post.merge();
OL.endif
Post.up(minimum.txt,1,5,Fields/Jpx, minimum.txt,1,7,Fields/Jpz);

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
