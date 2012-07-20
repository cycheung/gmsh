% comment lines start with "%" in any .ol file

% "radioButton" represent a boolean parameter or a flag
LOGFILES.radioButton(0,MetaModel/,''Output goes in .log files'');

% Flags to descibe model features that are activated or not
TENEUR.radioButton(1,Parameters/Model/1,"Account for variable water content"); 
CONVBC.radioButton(0,Parameters/Model/2,"Account for convection");

% Enumeration, i.e. a set of real values each associated with a label
SKINTYPE.number(1, Parameters/Model/4, ''Skin type''); 
SKINTYPE.addChoices(1,2);
SKINTYPE.addLabels(hairy, hairless);

% EPIDERMIS is determined by SKINTYPE
% Such dependency can be implemented with setValue
% The "setValue" statement overrules the value on the server. 
EPIDERMIS.number(0,Parameters/Model/5,''Epdidermis width [mm]'');
OL.if( OL.get(SKINTYPE) == 1)
EPIDERMIS.setValue(0.05);
OL.endif
OL.if( OL.get(SKINTYPE) == 2)
EPIDERMIS.setValue(0.12);
OL.endif

% onelab numbers
DERMIS.number(1.5,Parameters/Model/6,''Dermis width [mm]'');
WCONTENT.number(0.65,Parameters/Model/,''Water content []'');
BODYTEMP.number(305, Parameters/Model/,''Body temperature [K]'');
OVERTEMP.number(320, Parameters/Model/,''Thermal threshold fiber [K]'');

OL.if( OL.get(Parameters/Model/TENEUR) )
WCONTENT.setVisible(0);
OL.else
WCONTENT.setVisible(1);
OL.endif

% depending variables are defined with no value
% and this definition must then be completed by a "setValue" statement
ZSURF.number( , PostPro/);
ZSURF.setValue(OL.eval( (OL.get( DERMIS)+OL.get(EPIDERMIS))* 1e-3)); 

% "OL.get" return the value on server 
% of a parameter of type onelab::number or onelab::string
% "OL.eval" allows evaluating analytical expressions involving onelab::numbers

% The value of ZSURF is complemented with a list of choices
% which are the coordinates at which T will be monitored.
% The list of choice can be constructed element by element (as below) or by blocks: 
% param.addChoices(1,2,3); param.addChoices(7,12); 
% The 'value' of a parameter and the 'choices' can be evaluated independently
% according to the context and the needs.
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.001 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.049 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.100 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.150 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.200 * 1e-3) );

% Available LASER models, another enumeration
LASERTYPE.number(3, Parameters/Laser/1,''Laser type'');  
LASERTYPE.addChoices(1,2,3,4); 
LASERTYPE.addLabels(Applied temperature, Surface flux, Volume Flux, Controlled temperature);

APPLICTIME.number(0.05, Parameters/Laser/2, ''Application time [s]'');
ABSORPTION.number(2e4, Parameters/Laser/3, ''Absorption coefficient [1/m]'');
BEAMRADIUS.number(5, Parameters/Laser/4, ''Beam radius [mm]'');
REFLECTIVITY.number(0.0078, Parameters/Laser/5, ''Skin reflectivity []'');
LASERTEMP.number(360, Parameters/Laser/, ''Laser temperature [K]'');
LASERPOWER.number(15, Parameters/Laser/, ''Injected power [W]'');

% Visibility of the parameters in the onelab interactive window
% are controled with conditional statements
% so that only the relevant parameters appear.
OL.if( OL.get(LASERTYPE) == 1)
LASERTEMP.setVisible(1);
LASERPOWER.setVisible(0);
ABSORPTION.setVisible(0);
OL.endif
OL.if( OL.get(LASERTYPE) == 2)
LASERTEMP.setVisible(0);
LASERPOWER.setVisible(1);
ABSORPTION.setVisible(0);
OL.endif
OL.if( OL.get(LASERTYPE) == 3)
LASERTEMP.setVisible(0);
LASERPOWER.setVisible(1);
ABSORPTION.setVisible(1);
OL.endif
OL.if( OL.get(LASERTYPE) == 4)
LASERTEMP.setVisible(1);
LASERPOWER.setVisible(0);
ABSORPTION.setVisible(1);
OL.endif

% The metamodel is described as a list of clients in the "name.ol" file (this file)
% This metamodel has 6 clients

% syntax for clients
% OL.client name.Register([interf...|encaps...]{,cmdl{,wdir,{host{,rdir}}}}) ;

%-1)  Gmsh for meshing
Mesher.register(encapsulated);
Mesher.in( OL.get(Arguments/FileName).geo );
Mesher.args( OL.get(Arguments/FileName).geo);
Mesher.out( OL.get(Arguments/FileName).msh );
% Merge the mesh file if the metamodel is loaded by Gmsh
Mesher.merge( OL.get(Arguments/FileName).msh);

%-2) ElmerGrid converts the mesh for Elmer
ElmerGrid.register(interfaced);
ElmerGrid.in( OL.get(Arguments/FileName).msh);
ElmerGrid.args(14 2 OL.get(Arguments/FileName).msh -out mesh);
ElmerGrid.out( mesh/mesh.boundary );

%-3) ElmerSolver computes the thermal problem
Elmer.register(interfaced);
Elmer.in( ELMERSOLVER_STARTINFO.ol , OL.get(Arguments/FileName).sif.ol);
Elmer.out( solution.pos, temp.txt );

%-4) Post-processing with Gmsh and a script
Post.register(interfaced);
Post.in(solution.pos , script.opt.ol ); 
Post.args(solution.pos script.opt -);
Post.out(tempmin.txt, tempmax.txt, active*.txt, temp*.txt);
Post.up( tempmin.txt,-1,8,Solution/Tmin, tempmax.txt,-1,8,Solution/Tmax);

%-5) Display solution with a client Gmsh
Display.register(interfaced);
Display.in(solution.pos, script2.opt.ol, overheat.pos.opt.ol );
Display.out(overheat.pos );
Display.args( solution.pos script2.opt - );
Display.merge(overheat.pos);

%-6) Display solution curves with either gnuplot or matlab
POSTPRO.number(2, PostPro/,"Plot results with");
POSTPRO.addChoices(1,2); 
POSTPRO.addLabels(Matlab,Gnuplot);

Matlab.register(interfaced); 
Matlab.args(-nosplash -desktop -r plotMatlab);

Gnuplot.register(interfaced);
Gnuplot.in(temp.txt, plot.plt.ol);
Gnuplot.args(plot.plt );

OL.if( OL.get(POSTPRO) == 1)
Gnuplot.active(0);
Matlab.active(1);
OL.endif
OL.if( OL.get(POSTPRO) == 2)
Gnuplot.active(1);
Matlab.active(0);
OL.endif
