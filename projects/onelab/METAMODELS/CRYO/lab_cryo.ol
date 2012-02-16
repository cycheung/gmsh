%1) Client Gmsh pour faire le maillage
OL.client Mesher.Register(encapsulated);
OL.client Mesher/InputFiles.Set(OL.getValue(Arguments/FileName).geo);
%OL.client Mesher/LineOptions.Set(-v 0);

%2) Client ElmerGrid pour convertir le maillage pour Elmer
OL.client ElmerGrid.Register(interfaced);
OL.client ElmerGrid/PreLineOptions.Set(14 2);
OL.client ElmerGrid/InputFiles.Set(OL.getValue(Arguments/FileName).msh);
OL.client ElmerGrid/LineOptions.Set(-out mesh);
OL.client ElmerGrid/OutputFiles.Set(mesh/mesh.header);

%3) Client Elmer pour lancer la simulation avec Elmer
OL.client Elmer.Register(interfaced);
OL.client Elmer/InputFiles.Set(OL.getValue(Arguments/FileName).sif.ol);

%4)Client Postpro pour charger la solution calculée avec Elmer 
OL.client GmshMerge/InputFiles.Set(solution.pos);

%5)Client Postpro pour charger la solution calculée avec Elmer 
OL.client Post.Register(interfaced);
OL.client Post/InputFiles.Set(solution.pos,script.opt);
OL.client Post/LineOptions.Set("-");

%6)Client Postpro pour lancer Matlab/Gnuplot et dessiner la fonctionelle
OL.client Gnuplot.Register(interfaced);
OL.client Gnuplot/InputFiles.Set(plot.gnu);

%7)Client Postpro pour extraire la min et la max de la courbe
OL.client PostArray.List(fmax.txt,1,2,Solution/tmin, fmax.txt,1,8,Solution/fmin);


