MATERIAL.number(1, 3Material/0,"Name");
MATERIAL.valueLabels(1, "Steel", 2, "Wood", 3, "Alu", 4, "User defined");
YOUNG.number(, 3Material/, "Young modulus [GPa]");
POISSON.number(, 3Material/, "Poisson ratio []", 0:0.46:0.1);
POISSON.withinRange();
DENSITY.number(, 3Material/, "Mass density [kg/m3]");

OL.if(OL.get(MATERIAL) == 1)
  YOUNG.setValue(210);
  POISSON.setValue(0.3);
  DENSITY.setValue(7500);
OL.endif
OL.if(OL.get(MATERIAL) == 2)
  YOUNG.setValue(13);
  POISSON.setValue(0.45);
  DENSITY.setValue(900);
OL.endif
OL.if(OL.get(MATERIAL) == 3)
  YOUNG.setValue(70);
  POISSON.setValue(0.33);
  DENSITY.setValue(6000);
OL.endif
OL.if(OL.get(MATERIAL) == 4)
  YOUNG.setReadOnly(0);
  POISSON.setReadOnly(0);
  DENSITY.setReadOnly(0);
OL.endif

WEIGHT.radioButton(1, 3Material/0, "Account for weight");
LOOP.radioButton(0,,"Compute MT diagrams");

M.number(,9Results/2,"Moment - M [Nm]");
N.number(,9Results/3,"Traction - N [N]");
T.number(,9Results/4,"Shear - T [N]");
vmin.number(,9Results/,"min v(x) [m]"); 
vmax.number(,9Results/,"max v(x) [m]"); 
MAX.number(,9Results/1,"max |v(x)| [mm]");

vmin.setVisible(0);
vmax.setVisible(0);
M.setAttribute(Highlight, Coral);
N.setAttribute(Highlight, Coral);
T.setAttribute(Highlight, Coral);
MAX.setAttribute(Highlight, Coral);

X.number(0.10, 1Geometry/,"Cut location [m]");

OL.iftrue(LOOP)
X.setAttribute(Loop,3);
X.setAttribute(Graph,101000);
M.setAttribute(Graph,010000);
T.setAttribute(Graph,000100);
OL.else
  OL.if( OL.get(X, attrib.get(Loop)) == 3)
    X.setAttribute(Loop,0);
    X.setAttribute(Graph,0);
    M.setAttribute(Graph,0);
    T.setAttribute(Graph,0);
  OL.endif
OL.endif

OL.if( OL.get(X, choices.index()) == 0 )
  M.resetChoices();
  T.resetChoices();
OL.endif


#1) Client Gmsh pour le maillage initial
Mesher.register(native);
Mesher.in(OL.get(Arguments/FileName).geo);
Mesher.out(OL.get(Arguments/FileName).msh);
Mesher.run(OL.get(Arguments/FileName).geo -v 1);
Mesher.frontPage(OL.get(Arguments/FileName).geo,
                 OL.get(Arguments/FileName).msh);

OL.iftrue(LOOP)
X.range(0.10:OL.get(1Geometry/L)|15.00001);
OL.else
X.range(1e-4:OL.eval(OL.get(1Geometry/L)-1e-4):0.1);
X.withinRange();
OL.endif

#2) Client ElmerGrid pour convertir le maillage pour Elmer
ElmerGrid.register(interfaced);
ElmerGrid.in(OL.get(Arguments/FileName).msh);
ElmerGrid.out(mesh);
ElmerGrid.run(14 2 OL.get(Arguments/FileName).msh -out mesh);

#3) Client Elmer pour lancer la simulation avec Elmer
Elmer.register(encapsulated);
Elmer.in(ELMERSOLVER_STARTINFO.ol, beam.sif.ol);
Elmer.out( beam.pos, beam.txt);
Elmer.run();
Elmer.up(beam.txt,-1,1,9Results/vmin,
         beam.txt,-1,2,9Results/vmax);

#4) Post-processing with Gmsh and a script
Post.register(interfaced);
Post.in(beam.pos, MNT.pos.ol, script.opt.ol);
Post.out(MNT.txt);
Post.run(beam.pos MNT.pos -);
Post.up(MNT.txt,-1,8,9Results/2M,
        MNT.txt,-1,9,9Results/3N,
        MNT.txt,-1,10,9Results/4T);
Post.merge(RemoveViews.macro, beam.pos, script.opt);

MAX.setValue(OL.eval( 1e3*Max(abs(OL.get(vmin)), abs(OL.get(vmax))) ));

#OL.dump(aaa);