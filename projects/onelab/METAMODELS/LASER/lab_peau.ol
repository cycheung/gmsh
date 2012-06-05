LOGFILES.radioButton(0,MetaModel/,''Output goes in .log files'');

% Model parameters

TENEUR.radioButton(0,Parameters/Model/1,"Account for water content"); 
CONVBC.radioButton(0,Parameters/Model/2,"Account for convection");
BIOHEAT.radioButton(0,Parameters/Model/3,"Account for volume heat sources");

SKINTYPE.number(1, Parameters/Model/4, ''Skin type''); 
SKINTYPE.addChoices(1,2);
SKINTYPE.addLabels(hairy, hairless);

SKINWIDTH.number(0.05,Parameters/Model/5,''Skin width [mm]'');
//SKINWIDTH.setReadOnly(0);
OL.if( OL.get(SKINTYPE) == 1)
SKINWIDTH.setValue(0.05);
OL.endif
OL.if( OL.get(SKINTYPE) == 2)
SKINWIDTH.setValue(0.12);
OL.endif
DERMIS.number(1.5,Parameters/Model/6,''Dermis width [mm]'');

% z coordinates for post-processing curves
ZSURF0.number( OL.eval((OL.get(DERMIS)+OL.get(SKINWIDTH)-0.001)/1000), PostPro/);
ZSURF1.number( OL.eval((OL.get(DERMIS)+OL.get(SKINWIDTH)-0.050)/1000), PostPro/);
ZSURF2.number( OL.eval((OL.get(DERMIS)+OL.get(SKINWIDTH)-0.100)/1000), PostPro/);
ZSURF3.number( OL.eval((OL.get(DERMIS)+OL.get(SKINWIDTH)-0.150)/1000), PostPro/);
ZSURF4.number( OL.eval((OL.get(DERMIS)+OL.get(SKINWIDTH)-0.200)/1000), PostPro/);
ZSURF0.setReadOnly(1);
ZSURF1.setReadOnly(1);
ZSURF2.setReadOnly(1);
ZSURF3.setReadOnly(1);
ZSURF4.setReadOnly(1);


BEAMRADIUS.number(5, Parameters/Model/, ''Beam radius [mm]'');
WCONTENT.number(0.65,Parameters/Model/,''Water content []'');
BODYTEMP.number(310, Parameters/Model/, Body temperature [K]'');
OVERTEMP.number(320, Parameters/Model/, Maximum skin temperature [K]'');

%Available LASER models

LASERTYPE.number(1, Parameters/Laser/1,''Laser type'');  
LASERTYPE.addChoices(1,2,3); 
LASERTYPE.addLabels(Applied temperature, Surface flux, Volume Flux);

APPLICTIME.number(0.05, Parameters/Laser/, ''Application time [s]'');
LASERTEMP.number(360, Parameters/Laser/, ''Laser temperature [K]'');
LASERPOWER.number(15, Parameters/Laser/, ''Injected power [W]'');
ABSORPTION.number(2e4, Parameters/Laser/, ''Absorption coefficient [1/m]'');

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


% syntax for clients
% OL.client name.Register([interf...|encaps...]{,cmdl{,wdir,{host{,rdir}}}}) ;

%-1) Client Gmsh pour faire le maillage

Gmsh.register(encapsulated);
Gmsh.in( OL.get(Arguments/FileName).geo );
Gmsh.args( OL.get(Arguments/FileName).geo );
Gmsh.out( OL.get(Arguments/FileName).msh );
Gmsh.merge( OL.get(Arguments/FileName).msh);

%-2) Client ElmerGrid pour convertir le maillage pour Elmer
ElmerGrid.register(interfaced);
ElmerGrid.in( OL.get(Arguments/FileName).msh);
ElmerGrid.args(14 2 OL.get(Arguments/FileName).msh -out mesh);
ElmerGrid.out( mesh/mesh.boundary );

%-3) Client Elmer pour lancer la simulation avec Elmer
Elmer.register(interfaced);
Elmer.in( ELMERSOLVER_STARTINFO.ol , OL.get(Arguments/FileName).sif.ol);
Elmer.out( solution.pos, temp.txt );

%-4)Client Post pour lancer le script
Post.register(interfaced);
Post.in(solution.pos , script.opt.ol ); 
Post.args(solution.pos script.opt -);
Post.out(overheat.pos, temper.txt);
Post.merge(overheat.pos);

%GmshMerge.in(overheat.pos);

%-6) Client Postpro avec Matlab ou gnuplot

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


