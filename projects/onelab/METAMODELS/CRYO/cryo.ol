# global variables 

DISTANT.number(1,0Host/0,"host");
DISTANT.valueLabels(0,"local",1,"remote");
HOST.string(frhenrotte@arachnide_b,0Host/,"Remote login: user 'at' host"); 
REMDIR.string(ONELAB/CRYO/,0Host/,"Remote work directory");
ELMERDIR.string(,0Host/,"Remote bin Elmer directory");

RSYNCDELAY.number(1); RSYNCDELAY.setVisible(0);

HOST.setVisible(OL.get(DISTANT)); 
REMDIR.setVisible(OL.get(DISTANT));
ELMERDIR.setVisible(OL.get(DISTANT));

#1) Client Gmsh pour le maillage initial
Mesher.register(native);
Mesher.in( _OL.get(Arguments/FileName).geo);
Mesher.run(OL.get(Arguments/FileName).geo);
Mesher.out(OL.get(Arguments/FileName).msh);
Mesher.frontPage(OL.get(Arguments/FileName).geo, OL.get(Arguments/FileName).msh);

#2) Client ElmerGrid pour convertir le maillage pour Elmer
OL.iftrue(DISTANT)
  ElmerGrid.register(interfaced, OL.get(ELMERDIR)ElmerGrid, OL.get(HOST), OL.get(REMDIR));
OL.else
  ElmerGrid.register(interfaced);
OL.endif
ElmerGrid.in(_OL.get(Arguments/FileName).msh);
ElmerGrid.run(14 2 OL.get(Arguments/FileName).msh -out mesh);
ElmerGrid.out(mesh);

#3) Client Elmer pour lancer la simulation avec Elmer
OL.iftrue(DISTANT)
  Elmer.register(encapsulated, OL.get(ELMERDIR)ElmerSolver, OL.get(HOST),  OL.get(REMDIR));
OL.else
  Elmer.register(encapsulated);
OL.endif
Elmer.in(_ELMERSOLVER_STARTINFO.ol, _OL.get(Arguments/FileName).sif.ol);
Elmer.out(_solution.pos, _tempevol.txt);
Elmer.run();
Elmer.merge(solution.pos);

#4) Client Postpro pour executer un script gmsh
Post.register(interfaced);
Post.in(solution.pos, script.opt );
Post.run(solution.pos script.opt -);
Post.out(f.txt, fmax.txt);
Post.up(fmax.txt,-1,2,Solution/tmin, 
        fmax.txt,-1,8,Solution/fmin);

#6) Client Postpro pour extraire la min et la max de la courbe
#this client always run on the local machine
#Alternatively to gnuplot, matlab could be used. 
Gnuplot.register(interfaced);
Gnuplot.in(f.plt.ol, f.txt, fmax.txt);
Gnuplot.run(f.plt);


