OL.client Mesher.Register(encapsulated, gmsh);
OL.client Mesher/InputFiles.Set(OL.getValue(Arguments/FileName).geo.ol);
OL.client Mesher/LineOptions.Set(-v 0);

OL.client ElmerGrid.Register(interfaced,ElmerGrid 14 2);
OL.client ElmerGrid/InputFiles.Set(OL.getValue(Arguments/FileName).msh);
OL.client ElmerGrid/LineOptions.Set(-out meshdir);

OL.client Elmer.Register(interfaced);
OL.client Elmer/InputFiles.Set(OL.getValue(Arguments/FileName).sif.ol);

OL.client Post.Register(interfaced,gmsh);
OL.client Post/InputFiles.Set(solution.pos,script.opt);
OL.client Post/LineOptions.Set("-");

OL.client GmshMerge/InputFiles.Set(solution.pos);

OL.parameter TRANSIENT.number(1,Elmer/2,Transient simulation); TRANSIENT.AddChoices(0,1);
OL.iftrue(TRANSIENT)
	OL.client Gnuplot.Register(interfaced);
	OL.client Gnuplot/InputFiles.Set(f.plt);
OL.endif

OL.client PostArray.List(fmax.txt,1,2,tmin,fmax.txt,1,8,fmin);