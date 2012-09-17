# global variables 

DISTANT.radioButton(0,,"Flag remote host");

# DISTANT=0: local host
# DISTANT=1: remote host is localhost

WORKDIR.string();

OL.if(OL.get(DISTANT) == 0)
    HOST.string(); 
    BINDIR.string();
    REMDIR.string();
OL.endif
OL.if( OL.get(DISTANT) == 1)
# we choose localhost as the remote client.
    HOST.string(frhenrotte@localhost); 
    BINDIR.string(~/SCRIPTS/bash);
    REMDIR.string(/Users/frhenrotte/tmp/CRYO);
OL.endif

#1) Client Gmsh pour le maillage initial
OL.iftrue(DISTANT)
  Mesher.register(encapsulated, OL.get(BINDIR)/gmsh, OL.get(WORKDIR), OL.get(HOST), OL.get(REMDIR));
OL.else
  Mesher.register(encapsulated, gmsh);
OL.endif
Mesher.in( ./OL.get(Arguments/FileName).geo);
Mesher.args(OL.get(Arguments/FileName).geo);
Mesher.out(OL.get(Arguments/FileName).msh);

#2) Client ElmerGrid pour convertir le maillage pour Elmer
OL.iftrue(DISTANT)
  ElmerGrid.register(interfaced, OL.get(BINDIR)/ElmerGrid, OL.get(WORKDIR), OL.get(HOST), OL.get(REMDIR));
OL.else
  ElmerGrid.register(interfaced,ElmerGrid);
OL.endif
ElmerGrid.in( OL.get(Arguments/FileName).msh);
ElmerGrid.args(14 2 OL.get(Arguments/FileName).msh -out mesh);
ElmerGrid.out(mesh);

#3) Client Elmer pour lancer la simulation avec Elmer
OL.iftrue(DISTANT)
  Elmer.register(interfaced, OL.get(BINDIR)/ElmerSolver, OL.get(WORKDIR), OL.get(HOST),  OL.get(REMDIR));
OL.else
  Elmer.register(interfaced, ElmerSolver);
OL.endif
Elmer.in( ./ELMERSOLVER_STARTINFO.ol, ./OL.get(Arguments/FileName).sif.ol);
Elmer.out( solution.pos, tempevol.txt);
Elmer.merge(solution.pos);
# no args

#4) Client Postpro pour executer un script gmsh
OL.iftrue(DISTANT)
   Post.register(interfaced, OL.get(BINDIR)/ElmerSolver, OL.get(WORKDIR), OL.get(HOST),  OL.get(REMDIR));
OL.else
   Post.register(interfaced, gmsh);
OL.endif
Post.in( solution.pos , ./script.opt );
Post.args(solution.pos script.opt -);
Post.out( ./f.txt , ./fmax.txt);
Post.up(fmax.txt,-1,2,Solution/tmin, fmax.txt,-1,8,Solution/fmin);

#6) Client Postpro pour extraire la min et la max de la courbe
Gnuplot.register(interfaced, gnuplot);
Gnuplot.in(f.plt);
Gnuplot.args(f.plt);
