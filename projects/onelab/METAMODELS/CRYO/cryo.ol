# global variables 

DISTANT.number(0, 0Host/0, "host");
DISTANT.valueLabels(0,"local",
                    1,"remote");
OL.iftrue(DISTANT)
#ElmerGrid.remote(frhenrotte@arachnide_b, ONELAB/CRYO/);
#Elmer.remote(frhenrotte@arachnide_b, ONELAB/CRYO/);
ElmerGrid.remote();
Elmer.remote();
OL.endif

#1) Client Gmsh pour le maillage initial
Mesher.register(native);
Mesher.in( _OL.get(Arguments/FileName).geo);
Mesher.run(OL.get(Arguments/FileName).geo);
Mesher.out(OL.get(Arguments/FileName).msh);
Mesher.frontPage(OL.get(Arguments/FileName).geo, OL.get(Arguments/FileName).msh);

#2) Client ElmerGrid pour convertir le maillage pour Elmer
ElmerGrid.register(interfaced);
ElmerGrid.in( _OL.get(Arguments/FileName).msh);
ElmerGrid.run(14 2 OL.get(Arguments/FileName).msh -out mesh);
ElmerGrid.out(mesh);

#3) Client Elmer pour lancer la simulation avec Elmer
Elmer.register(encapsulated);
Elmer.in( _ELMERSOLVER_STARTINFO.ol, _OL.get(Arguments/FileName).sif.ol);
Elmer.out( _solution.pos, _tempevol.txt);
Elmer.run();
Elmer.merge(solution.pos);

tmin.number(,Solution/,"Time when f(t) is minimum");
fmin.number(,Solution/,"Minimum of f(t)");

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



