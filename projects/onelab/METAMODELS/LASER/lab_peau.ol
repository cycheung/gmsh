# Onelab commands start with the tag "OL."
# Comment lines start with the tag "#"
# These tags are the default values.
# They can be modified to accomodate various client syntaxes.
# This is done with the sentence e.g.: onelab.tags(/,//);
# Defaults are restored with
# onelab.tags(); onelab.tags(,); or onelab.tags(OL.,#);


# The metamodel is described as a list of clients 
# in the "name.ol" file (this file)
# In this case, the metamodel has 6 clients

#-1)  Gmsh for meshing
Mesher.register(native);
Mesher.in( OL.get(Arguments/FileName).geo );
Mesher.run( OL.get(Arguments/FileName).geo );
Mesher.out( OL.get(Arguments/FileName).msh );
# Tell ONELAB to show geometry + mesh as initial view
# Mesher.preCompute();
# Mesher.frontPage(OL.get(Arguments/FileName).geo, OL.get(Arguments/FileName).msh);
# Tell ONELAB to show geometry as initial view
Mesher.frontPage(OL.get(Arguments/FileName).geo);

# Enumeration, i.e. a set of real values each associated with a label
SKINTYPE.number(2, Parameters/Skin/1,"Skin type"); 
SKINTYPE.valueLabels(1,"hairy (thin)", 
                     2,"glabrous (thick)");
# Numbers in pathes allow to sort parameters in the onelab window.
# SKINTYPE will be the 1 paraemeter in the subtree /Parameters/Skin

# The thickness of epidermis (parameter EPIDERMIS) is determined 
# by the value of SKINTYPE, i.e. it is a depending variable.
# In .ol files, depending variables are declared with no value 
# (the value slot is left empty)
# and the incomplete declaration is completed by a "setValue" statement
# In this case, EPIDERMIS was defined in lab_peau.geo

OL.if( OL.get(SKINTYPE) == 1)
Parameters/Skin/EPIDERMIS.setValue(0.05);
OL.endif
OL.if( OL.get(SKINTYPE) == 2)
Parameters/Skin/EPIDERMIS.setValue(0.12);
OL.endif
# The "setValue" statement overrules the value on server.

# other parameters of the model
WCONTENT.number(0.65,Parameters/Skin/,"Water content []");
BODYTEMP.number(37, Parameters/Skin/,"Body temperature [C]");
OVERTEMP.number(47, Parameters/Skin/,"Thermal threshold fiber [C]");
REFLECTIVITY.number(0.0078, Parameters/Skin/, "Skin reflectivity []");

# Flags to describe model features that are activated or not
TENEUR.radioButton(1,Parameters/Skin/3,"Account for variable water content");
OL.iftrue(TENEUR)
WCONTENT.setVisible(0);
OL.else
WCONTENT.setVisible(1);
OL.endif

CONVBC.radioButton(0,Parameters/Skin/4,"Account for convection");
HCONV.number(100, Parameters/Skin/5, "Convection coefficient [W/(Km)]");
TAMBIANT.number(20, Parameters/Skin/6, "Ambiant temperature [C]");
OL.iftrue(CONVBC) 
HCONV.setVisible(1);
OL.else
HCONV.setVisible(0);
OL.endif

# Available LASER models, another enumeration
LASERTYPE.number(2, Parameters/Laser/1,"Laser type");  
LASERTYPE.valueLabels(
    1,"Imposed temperature", 
    2,"Imposed flux", 
    3,"Controlled temperature");

LASERSHAPE.number(1, Parameters/Laser/2,"Laser shape");  
LASERSHAPE.valueLabels(1,"Gaussian", 2,"Flat-top");

# Parameters describing the laser stimulator
STIMTIME.number(0.110, Parameters/Laser/, "Stimulus duration [s]");
ABSORPTION.number(2e4, Parameters/Laser/, "Absorption coefficient [1/m]");
LASERTEMP.number(50, Parameters/Laser/, "Imposed temperature [C]");
LASERPOWER.number(4.0, Parameters/Laser/, "Power [W]");
TSKINFILE.string(,Parameters/Laser/, "Imposed temp. file");
QSKINFILE.string(,Parameters/Laser/, "Imposed laser power file");
TSKINFILE.addChoices(tskin.sif);
QSKINFILE.addChoices(pin.sif);

# Visibility of the parameters in the onelab interactive window
# can be controlled with conditional statements
# so that only the relevant parameters appear.
OL.if( OL.get(LASERTYPE) == 1)
LASERTEMP.setVisible(1);
TSKINFILE.setVisible(1);
LASERPOWER.setVisible(0);
QSKINFILE.setVisible(0);
ABSORPTION.setVisible(0);
LASERSHAPE.setVisible(0);
REFLECTIVITY.setVisible(0);
OL.endif
OL.if( OL.get(LASERTYPE) == 2)
LASERTEMP.setVisible(0);
TSKINFILE.setVisible(0);
LASERPOWER.setVisible(1);
QSKINFILE.setVisible(1);
ABSORPTION.setVisible(1);
LASERSHAPE.setVisible(1);
REFLECTIVITY.setVisible(1);
OL.endif
OL.if( OL.get(LASERTYPE) == 3)
LASERTEMP.setVisible(1);
TSKINFILE.setVisible(0);
LASERPOWER.setVisible(1);
QSKINFILE.setVisible(0);
ABSORPTION.setVisible(1);
LASERSHAPE.setVisible(1);
REFLECTIVITY.setVisible(1);
OL.endif

PROBETIME.number(OL.get(STIMTIME), PostPro/, "Probe time [s]");

ZSURF.number( , PostPro/,"Z coordinates");
ZSURF.setValue(OL.eval( (OL.get(Parameters/Skin/DERMIS)+OL.get(Parameters/Skin/EPIDERMIS))* 1e-3)); 

# "OL.get" return the value on server 
# of a parameter of type onelab::number or onelab::string
# "OL.eval" allows evaluating analytical expressions involving onelab::numbers

# The parameter ZSURF is attributed a list of choices
# which are the coordinates at which the temperature will be monitored.
# The list of choice can be constructed element by element (as below) 
# or by blocks: param.addChoices(1,2,3); param.addChoices(7,12); 
# The 'value' of a parameter and the 'choices' are different concepts
# used independently according to the context.
ZSURF.resetChoices();
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.0001 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.0125 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.0250 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.0375 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.0501 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.0625 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.0750 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.0875 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.1000 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.1125 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.1250 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.1375 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.1500 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.1625 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.1750 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.1875 * 1e-3) );
ZSURF.addChoices( OL.eval( OL.get(ZSURF) - 0.2000 * 1e-3) );

VOLUME.number(,Solution/,"Vol. above threshold [mm3]");

#-2) ElmerGrid converts the mesh for Elmer
ElmerGrid.register(interfaced);
ElmerGrid.in( OL.get(Arguments/FileName).msh);
ElmerGrid.out( mesh/mesh.boundary );
ElmerGrid.run(14 2 OL.get(Arguments/FileName).msh -out mesh);


#-3) ElmerSolver computes the thermal problem
Elmer.register(encapsulated);
Elmer.in( ELMERSOLVER_STARTINFO.ol, OL.get(Arguments/FileName).sif.ol);
Elmer.out( solution.pos, temp.txt );
Elmer.run( );
Elmer.merge(solution.pos);


#-4) Post-processing with Gmsh and a script
Post1.register(interfaced);
Post1.in(solution.pos , script.opt.ol );
Post1.out(temp0.txt, activeMax.txt);
Post1.run(solution.pos script.opt -);

#-5) Further post-processing with Gmsh and a script
Post2.register(interfaced);
Post2.in(solution.pos, script2.opt.ol);
Post2.out(volume.txt, aboveThres.pos );
Post2.run( solution.pos script2.opt - );
Post2.merge(aboveThres.pos);
Post2.up( volmax.txt,-1,8,Solution/VOLUME);

#-6) Display solution curves with either gnuplot or matlab
POSTPRO.number(2, PostPro/,"Plot results with");
POSTPRO.valueLabels(1, "Matlab",2, "Gnuplot");
POSTPRO.setVisible(0);

OL.if( OL.get(POSTPRO) == 1)
Matlab.register(interfaced); 
Matlab.run(-nosplash -desktop -r plotMatlab);
OL.endif
OL.if( OL.get(POSTPRO) == 2)
Gnuplot.register(interfaced);
Gnuplot.in(temp.txt, plot.plt.ol);
Gnuplot.run(plot.plt );
OL.endif

# Dump the ONELAB database in a file named zzz
# OL.dump(zzz);