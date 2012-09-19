# global variables 

DISTANT.number(1,0Host/0,"host");
DISTANT.valueLabels(0,"local",1,"remote");

HOST.string(,0Host/,"Remote login: user 'at' host"); HOST.layout();
REMDIR.string(,0Host/,"Remote work directory"); REMDIR.layout();
ELMERDIR.string(,0Host/,"Remote Elmer directory"); ELMERDIR.layout();

HOST.setVisible(OL.get(DISTANT)); 
REMDIR.setVisible(OL.get(DISTANT));
ELMERDIR.setVisible(OL.get(DISTANT));

#1) Client Gmsh pour le maillage initial
Mesher.register(encapsulated, gmsh);
Mesher.in( ./OL.get(Arguments/FileName).geo);
Mesher.args(OL.get(Arguments/FileName).geo);
Mesher.out(OL.get(Arguments/FileName).msh);

#2) Client ElmerGrid pour convertir le maillage pour Elmer
OL.iftrue(DISTANT)
  ElmerGrid.register(interfaced, OL.get(ELMERDIR)ElmerGrid, OL.get(HOST), OL.get(REMDIR));
OL.else
  ElmerGrid.register(interfaced, ElmerGrid);
OL.endif
ElmerGrid.in(./OL.get(Arguments/FileName).msh);
ElmerGrid.args(14 2 OL.get(Arguments/FileName).msh -out mesh);
ElmerGrid.out(mesh);

#3) Client Elmer pour lancer la simulation avec Elmer
OL.iftrue(DISTANT)
  Elmer.register(interfaced, OL.get(ELMERDIR)ElmerSolver, OL.get(HOST),  OL.get(REMDIR));
OL.else
  Elmer.register(interfaced, ElmerSolver);
OL.endif
Elmer.in(./ELMERSOLVER_STARTINFO.ol, ./OL.get(Arguments/FileName).sif.ol);
Elmer.out( ./solution.pos, ./tempevol.txt);
Elmer.merge(solution.pos);
# no args

#4) Client Postpro pour executer un script gmsh
Post.register(interfaced, gmsh);
Post.in( solution.pos, script.opt );
Post.args(solution.pos script.opt -);
Post.out(f.txt, fmax.txt);
Post.up(fmax.txt,-1,2,Solution/tmin, 
        fmax.txt,-1,8,Solution/fmin);

#6) Client Postpro pour extraire la min et la max de la courbe
#this client always run on the local machine
#Alternatively to gnuplot, matlab could be used. 
Gnuplot.register(interfaced, gnuplot);
Gnuplot.in(f.plt.ol, f.txt, fmax.txt);
Gnuplot.args(f.plt);


