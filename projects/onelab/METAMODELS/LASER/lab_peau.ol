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

TAGSIMU.string(,,"Simulation tag"); 
TAGSIMU.setReadOnly(0); 
RESTORE.radioButton(0,Solution/0,"Restore");

#-1)  Gmsh for meshing
Mesher.register(native);
Mesher.in( OL.get(Arguments/FileName).geo );
Mesher.run( OL.get(Arguments/FileName).geo);
Mesher.out( OL.get(Arguments/FileName).msh );
# merge .geo file as initial view
Mesher.frontPage(OL.get(Arguments/FileName).geo);

# Enumeration, i.e. a set of real values each associated with a label
SKINTYPE.number(1, Parameters/Skin/1,"Skin type"); 
SKINTYPE.valueLabels(1,"hairy (thin)", 
                     2,"glabrous (thick)");
# Numbers in pathes allow to sort parameters in the onelab window.
# SKINTYPE will be the 1st parameter in the subtree /Parameters/Skin

# The thickness of epidermis (parameter EPIDERMIS) is determined 
# by the value of SKINTYPE, i.e. it is a depending variable.
# In .ol files, depending variables are declared with no value 
# (the value slot is left empty)
# and the incomplete declaration is completed by a "setValue" statement
# In this case, EPIDERMIS was defined in lab_peau.geo

OL.if( OL.get(SKINTYPE) == 1)
 Parameters/Skin/2EPIDERMIS.setValue(0.05);
OL.endif
OL.if( OL.get(SKINTYPE) == 2)
 Parameters/Skin/2EPIDERMIS.setValue(0.12);
OL.endif

# The "setValue" statement overrules the value on server.

# other parameters of the model
REFLECTIVITY.number(0.0078, Parameters/Skin/4, "Skin reflectivity []");
TCORE.number(37, Parameters/Skin/5,"Core temperature [C]");
TBASE.number(30, Parameters/Skin/6, "Baseline temperature [C]");

# Flags describing model features that are activated or not
TENEUR.radioButton(1,Parameters/Skin/8,"Account for variable water content");
CONVBC.radioButton(0,Parameters/Skin/9,"Account for convection");

WCONTENT.number(0.65,Parameters/Skin/,"Water content []");
OL.iftrue(TENEUR)
WCONTENT.setVisible(0);
OL.else
WCONTENT.setVisible(1);
OL.endif

HCONV.number(100, Parameters/Skin/, "Convection coefficient [W/(Km)]");
TAMBIANT.number(20, Parameters/Skin/, "Ambiant temperature [C]");
OL.iftrue(CONVBC) 
HCONV.setVisible(1);
TAMBIANT.setVisible(1);
OL.else
HCONV.setVisible(0);
TAMBIANT.setVisible(0);
OL.endif

# Available LASER models, another enumeration
LASERSHAPE.number(2, Parameters/Laser/1,"Laser shape");  
LASERSHAPE.valueLabels(1,"Gaussian",	    
                       2,"Flat-top");

# Parameters describing the laser stimulator
ABSORPTION.number(2e4, Parameters/Laser/3, "Absorption coefficient [1/m]");


LASERSIMU.number(3, Parameters/Laser/4,"Laser type");  
LASERSIMU.valueLabels(
    1,"Imposed power", 
    2,"Imposed power file", 
    3,"Controlled temperature",
    4,"Controlled temperature file");
LASERTYPE.number(1, Parameters/Laser/,"Laser category");
OL.if( OL.get(LASERSIMU) == 1)
Parameters/Laser/LASERTYPE.setValue(1);
OL.endif
OL.if( OL.get(LASERSIMU) == 2)
Parameters/Laser/LASERTYPE.setValue(1);
OL.endif
OL.if( OL.get(LASERSIMU) == 3)
Parameters/Laser/LASERTYPE.setValue(2);
OL.endif
OL.if( OL.get(LASERSIMU) == 4)
Parameters/Laser/LASERTYPE.setValue(2);
OL.endif
LASERTYPE.setVisible(0);

LASERPOWER.number(10.0, Parameters/Laser/5, "Power [W]");
LASERTEMP.number(50, Parameters/Laser/5, "Target temperature T_hot [C]");
STIMTIME.number(0.10, Parameters/Laser/6, "Pulse Duration [s]");

TEMP1.number(30, Parameters/Laser/7, "Conditioning temperature T_cond [C]");
TIME1.number(0.0, Parameters/Laser/8, "Conditioning time [s]");
RAMPTIME.number(0.01, Parameters/Laser/9, "Ramp time [s]");

TimeStep.number(1, Parameters/Elmer/2,"Time step [ms]");
TimeEnd.number(0.25,Parameters/Elmer/3,"Simulation end time [s]");

TSKINFILE.string(,Parameters/Laser/6, "Imposed temp. file");
QSKINFILE.string(,Parameters/Laser/6, "Imposed laser power file");
TSKINFILE.addChoices(Timposed.sif, tskin.sif, TMoureaux.sif);
QSKINFILE.addChoices(pin.sif);

# Visibility of the parameters in the onelab interactive window
# can be controlled with conditional statements
# so that only the relevant parameters appear.
OL.if( OL.get(LASERSIMU) == 1)
LASERTEMP.setVisible(0);
TSKINFILE.setVisible(0);
LASERPOWER.setVisible(1);
QSKINFILE.setVisible(0);
TSKINFILE.setVisible(0);
STIMTIME.setVisible(1);
TEMP1.setVisible(0);
TIME1.setVisible(0);
RAMPTIME.setVisible(0);
OL.endif

OL.if( OL.get(LASERSIMU) == 2)
LASERTEMP.setVisible(0);
TSKINFILE.setVisible(0);
LASERPOWER.setVisible(0);
QSKINFILE.setVisible(1);
TSKINFILE.setVisible(0);
STIMTIME.setVisible(0);
TEMP1.setVisible(0);
TIME1.setVisible(0);
RAMPTIME.setVisible(0);
OL.endif

OL.if( OL.get(LASERSIMU) == 3)
LASERTEMP.setVisible(1);
TSKINFILE.setVisible(0);
LASERPOWER.setVisible(1);
QSKINFILE.setVisible(0);
TSKINFILE.setVisible(0);
STIMTIME.setVisible(1);
TEMP1.setVisible(1);
TIME1.setVisible(1);
RAMPTIME.setVisible(1);
OL.endif

OL.if( OL.get(LASERSIMU) == 4)
LASERTEMP.setVisible(0);
TSKINFILE.setVisible(1);
LASERPOWER.setVisible(1);
QSKINFILE.setVisible(0);
TSKINFILE.setVisible(1);
STIMTIME.setVisible(0);
TEMP1.setVisible(0);
TIME1.setVisible(0);
RAMPTIME.setVisible(0);
OL.endif

ABSORPTION.setVisible(1);
LASERSHAPE.setVisible(1);
REFLECTIVITY.setVisible(1);

# The "Probe time" is the time at which temperature distributions
# and activation characteristics (volume and area) are evaluated.
# ONELAB can not always guess this important parameter 

PULSEEND.number(, PostPro/1, "Guessed pulse endtime [s]");
OL.if( OL.get(LASERSIMU) == 1)
  PULSEEND.setValue(OL.get(STIMTIME));
OL.endif
OL.if(OL.get(LASERSIMU) == 2) 
  PULSEEND.setValue(0);
OL.endif
OL.if( OL.get(LASERSIMU) == 3)
  PULSEEND.setValue(OL.eval(OL.get(TIME1)+OL.get(RAMPTIME)+OL.get(STIMTIME)));
OL.endif
OL.if(OL.get(LASERSIMU) == 4)
  PULSEEND.setValue(0);
OL.endif

PROBETIME.number(OL.get(PULSEEND), PostPro/2, "Probe time [s]");
OL.if( OL.get(PROBETIME) > OL.get(TimeEnd) )
  PROBETIME.setValue( OL.eval(OL.get(TimeEnd)*0.99) );
  PROBETIME.setReadOnly(0);
OL.endif

OVERTEMP.number(39.5, PostPro/5,"Thermal threshold fiber [C]");

SKINWIDTH.number(, PostPro/3,"Skin width [mm]");
PostPro/3SKINWIDTH.setValue( OL.eval(OL.get(Parameters/Skin/3DERMIS) + OL.get(Parameters/Skin/2EPIDERMIS)));

# The parameter ZSURF is attributed a list of choices
# which are the coordinates at which the temperature will be monitored
# in the top right diagram (pdf output).

ZSURF.number(0, PostPro/4,"Probe depth [mm]") ; #just below skin surface
ZSURF.resetChoices();
ZSURF.addChoices( OL.eval( (OL.get(PostPro/3SKINWIDTH) - 1e-4)));
ZSURF.addChoices( OL.eval( OL.get(Parameters/Skin/3DERMIS) - 1e-4));
ZSURF.addChoices( OL.eval( OL.get(PostPro/3SKINWIDTH) - OL.get(PostPro/4ZSURF) - 1e-04));

# Highlight post-processing parmeters that need be set by the user
PROBETIME.setAttribute(Highlight,AliceBlue);
OVERTEMP.setAttribute(Highlight,AliceBlue);
ZSURF.setAttribute(Highlight,AliceBlue);

# "OL.get" return the value on server 
#  of a parameter of type onelab::number or onelab::string
# "OL.eval" allows evaluating analytical expressions involving onelab::numbers

# The list of choice of a parameter can be constructed element by element
# or by blocks: param.addChoices(1,2,3); param.addChoices(7,12); 
# The 'value' of a parameter and the 'choices' are different concepts
# used independently according to the context.

MAXTEMPTIME.number(,Solution/,"Max temp. time at interface [s]");
TEMPINTERF.number(,Solution/,"Max temperature at interface [C]");
UADELTA.number(,Solution/,"Area A.T. at interface [mm2]");
VOLUME.number(,Solution/,"Skin volume A.T. [mm3]");
VOLDERM.number(,Solution/,"Derm volume A.T. [mm3]");
OL.if( OL.get(Metamodel/Action) == 1 ) # "check"
  VOLUME.resetChoices();
  VOLDERM.resetChoices();
  UADELTA.resetChoices();
OL.endif

#-2) ElmerGrid converts the mesh for Elmer
ElmerGrid.register(interfaced);
OL.ifntrue(Solution/0RESTORE)
  ElmerGrid.in( OL.get(Arguments/FileName).msh);
  ElmerGrid.out( mesh/mesh.boundary );
  ElmerGrid.run(14 2 OL.get(Arguments/FileName).msh -out mesh);
OL.endif

#-3) ElmerSolver computes the thermal problem
Elmer.register(encapsulated);
OL.ifntrue(Solution/0RESTORE)
  Elmer.in( ELMERSOLVER_STARTINFO.ol, OL.get(Arguments/FileName).sif.ol);
  Elmer.out( solutionOL.get(TAGSIMU).pos, tempOL.get(TAGSIMU).txt );
  Elmer.run( );
OL.endif

#-4a) Post-processing with Gmsh and a script
Post0.register(interfaced);
Post0.in(solutionOL.get(TAGSIMU).pos , script0.opt.ol );
Post0.out(tempMaxInterface.txt);
Post0.run(solutionOL.get(TAGSIMU).pos script0.opt -);
Post0.up(tempMaxInterface.txt,-1,2,Solution/MAXTEMPTIME,
         tempMaxInterface.txt,-1,8,Solution/TEMPINTERF,
         Adelta.txt,-1,8,Solution/UADELTA);

#-4b) Post-processing with Gmsh and a script
Post1.register(interfaced);
Post1.in(solutionOL.get(TAGSIMU).pos , script.opt.ol );
Post1.out(temp0.txt, activeMax.txt);
Post1.run(solutionOL.get(TAGSIMU).pos script.opt -);

#-5) Further post-processing with Gmsh and a script
Post2.register(interfaced);
Post2.in(solutionOL.get(TAGSIMU).pos, script2.opt.ol);
Post2.out(volume.txt, aboveThres.pos );
Post2.run( solutionOL.get(TAGSIMU).pos script2.opt - );
Post2.up(volmax.txt,-1,8,Solution/VOLUME,
         voldermmax.txt,-1,8,Solution/VOLDERM);

MERGE.radioButton(1,PostPro/,"Color plot display");
OL.iftrue(MERGE)
  Post2.merge(solutionOL.get(TAGSIMU).pos, aboveThres.pos);
OL.endif

#-6) Display solution curves with either gnuplot or matlab

Gnuplot.register(interfaced);
Gnuplot.in(tempOL.get(TAGSIMU).txt, volume.txt, volderm.txt, activeMax.txt, duration.txt,
           plot.plt.ol, 
           plot1.plt.ol, 
           plot2.plt.ol, 
           plot3.plt.ol, 
           plot4.plt.ol);
Gnuplot.run(plot.plt);

# Dump the ONELAB database in a file named zzz
#OL.dump(zzz);
#OL.show(Gmsh/MergedGeo, Arguments/WorkingDir);
OL.show(TAGSIMU);

OL.clientStatus();

