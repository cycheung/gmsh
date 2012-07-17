

DISPLAY.radioButton(0,Options/,"Display functional at each step");

Hx.number(0.5, Fields/,"Hx"); 
Hx.addChoices( 0, 0.25, 0.5, 0.75, 1, 2, 1.01,
               0.76, 0.51, 0.26, 0.01, -0.25, -0.5, -1,
	      -0.51, 0.02);

%Hx.addChoices( 0, 0.15, 0.5, 1, 1.5, 2,
%               1.51, 1.01, 0.51, 0.01, -0.5, -1, -1.5,
%	      -1.01, -0.51, -0.01, 0.49, 0.26, 0.02);

Hz.number(0, Fields/,"Hz");
Jpx.number(0, Fields/,"Jpx");
Jpz.number(0, Fields/,"Jpz");

% Client 1: 

Post.register (interfaced, gmsh);
Post.in ( circle.msh , script.ol );

Post.out( minimum.txt );
Post.args(circle.msh script -);
OL.iftrue(DISPLAY)
  Post.merge(Functional.pos);
OL.endif
Post.up(minimum.txt,1,5,Fields/Jpx, minimum.txt,1,7,Fields/Jpz);

% Client 2:

FIRSTTIME.radioButton(1); FIRSTTIME.setVisible(0);

Shell.register(interfaced, sh);
Shell.in( results.sh.ol );
Shell.args( results.sh ); 
%Shell.register(interfaced, echo);
%Shell.args("OL.get(Fields/Hx) OL.get(Fields/Jpx)" >> results.txt);
OL.if( OL.get(Fields/Hx) == 0 )
   Shell.out(results.txt);
OL.else
   Shell.out();
OL.endif

% Client 3:

OL.if( OL.get(Fields/Hx) == 0.02)
  Gnuplot.register(interfaced,gnuplot);
  Gnuplot.in(results.txt, plot.plt);
  Gnuplot.args(plot.plt);
OL.endif
