% comment lines start with "%" in any .ol file

% "radioButton" represent a boolean parameter or a flag
LOGFILES.radioButton(0,MetaModel/,''Output goes in .log files'');

% Flags to descibe model features that are activated or not
TENEUR.radioButton(0,Parameters/Model/1,"Account for water content"); 
CONVBC.radioButton(0,Parameters/Model/2,"Account for convection");
BIOHEAT.radioButton(0,Parameters/Model/3,"Account for volume heat sources");

% Enumeration, i.e. a set of real values each associated with a label
SKINTYPE.number(1, Parameters/Model/4, ''Skin type''); 
SKINTYPE.addChoices(1,2);
SKINTYPE.addLabels(hairy, hairless);

% SKINWIDTH is determined by SKINTYPE
% Such dependency can be implemented with setValue
% The "setValue" sentence overrules the value on the server. 
SKINWIDTH.number(0.05,Parameters/Model/5,''Skin width [mm]'');
OL.if( OL.get(SKINTYPE) == 1)
SKINWIDTH.setValue(0.05);
OL.endif
OL.if( OL.get(SKINTYPE) == 2)
SKINWIDTH.setValue(0.12);
OL.endif

% onelab numbers
DERMIS.number(1.5,Parameters/Model/6,''Dermis width [mm]'');
BEAMRADIUS.number(5, Parameters/Model/, ''Beam radius [mm]'');
WCONTENT.number(0.65,Parameters/Model/,''Water content []'');
BODYTEMP.number(310, Parameters/Model/, Body temperature [K]'');
OVERTEMP.number(320, Parameters/Model/, Maximum skin temperature [K]'');

% z coordinates for post-processing curves
ZSURF0.number( OL.eval((OL.get(DERMIS)+OL.get(SKINWIDTH)-0.001)/1000), PostPro/);
ZSURF1.number( OL.eval((OL.get(DERMIS)+OL.get(SKINWIDTH)-0.050)/1000), PostPro/);
ZSURF2.number( OL.eval((OL.get(DERMIS)+OL.get(SKINWIDTH)-0.100)/1000), PostPro/);
ZSURF3.number( OL.eval((OL.get(DERMIS)+OL.get(SKINWIDTH)-0.150)/1000), PostPro/);
ZSURF4.number( OL.eval((OL.get(DERMIS)+OL.get(SKINWIDTH)-0.200)/1000), PostPro/);
% "OL.get" return the value on server of a parameter of type onelab::number or onelab::string
% "OL.eval" allows evaluating analytical expressions involving onelab::numbers

% In the above definition, ZSURFx is determined by DERMIS and SKINWIDTH
% With the readOnly flag set to true, 
% parsing the defnition ZSURFx.number(...) will always overrules the value on server
% whereas with the readOnly flag set to false, the value on server would be retained.
ZSURF0.setReadOnly(1); 
ZSURF1.setReadOnly(1);
ZSURF2.setReadOnly(1);
ZSURF3.setReadOnly(1);
ZSURF4.setReadOnly(1);

% The lines:
% ZSURF0.number( OL.eval((OL.get(DERMIS)+OL.get(SKINWIDTH)-0.001)/1000), PostPro/);
% ZSURF0.setReadOnly(1);
% are equivalent to:
% ZSURF0.number( 0, PostPro/);
% ZSURF0.setValue( OL.eval((OL.get(DERMIS)+OL.get(SKINWIDTH)-0.001)/1000));
% It is a matter of choice to use "setValue" or "setReadOnly"

% Available LASER models, another enumeration
LASERTYPE.number(1, Parameters/Laser/1,''Laser type'');  
LASERTYPE.addChoices(1,2,3); 
LASERTYPE.addLabels(Applied temperature, Surface flux, Volume Flux);

APPLICTIME.number(0.05, Parameters/Laser/, ''Application time [s]'');
LASERTEMP.number(360, Parameters/Laser/, ''Laser temperature [K]'');
LASERPOWER.number(15, Parameters/Laser/, ''Injected power [W]'');
ABSORPTION.number(2e4, Parameters/Laser/, ''Absorption coefficient [1/m]'');

% Visibility of the parameters in the onelab interactive window
% are controled with conditional statements
% so that only the relevant parameters appear.
OL.if( OL.get(LASERTYPE) == 1)
LASERTEMP.setVisible(1);
LASERPOWER.setVisible(0);
OL.endif
OL.if( OL.get(LASERTYPE) == 2)
LASERTEMP.setVisible(0);
LASERPOWER.setVisible(1);
OL.endif
OL.if( OL.get(LASERTYPE) == 3)
LASERTEMP.setVisible(0);
LASERPOWER.setVisible(1);
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
Post.out(overheat.pos, temper.txt);
Post.merge(overheat.pos);

%-5) Display solution curves with either gnuplot or matlab
POSTPRO.number(2, PostPro/,"Plot results with");
POSTPRO.addChoices(1,2); 
POSTPRO.addLabels(Matlab,Gnuplot);
OL.if( OL.get(POSTPRO) == 1)
Matlab.register(interfaced); 
Matlab.args(-nosplash -desktop -r plotMatlab);
OL.endif
OL.if( OL.get(POSTPRO) == 2)
Gnuplot.register(interfaced);
Gnuplot.args(plot.plt );
OL.endif

%-6) Display solution with a client Gmsh if the metamodel runs in standalone
% i.e. not as called from Gmsh
Display.register(interfaced);
Display.args(OL.get(Arguments/FileName).msh overheat.pos);
OL.iftrue(HasGmsh)
Display.active(0);
OL.endif