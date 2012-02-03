OL.client Mesher.Register(encapsulated);
OL.client Mesher/InputFiles.Set(OL.getValue(Arguments/FileName).geo);
OL.client Mesher/LineOptions.Set(-v 0);

OL.client ElmerGrid.Register(interfaced);
OL.client ElmerGrid/PreLineOptions.Set(14 2);
OL.client ElmerGrid/InputFiles.Set(OL.getValue(Arguments/FileName).msh);
OL.client ElmerGrid/LineOptions.Set(-out meshdir);
OL.client ElmerGrid/OutputFiles.Set(meshdir/mesh.header);

OL.client Elmer.Register(interfaced);
OL.client Elmer/InputFiles.Set(OL.getValue(Arguments/FileName).sif.ol);

OL.client Post.Register(interfaced);
OL.client Post/InputFiles.Set(solution.pos,script.opt);
OL.client Post/LineOptions.Set("-");

OL.client GmshMerge/InputFiles.Set(solution.pos);

OL.parameter TRANSIENT.number(0,,Transient simulation); TRANSIENT.AddChoices(0,1);

OL.client Gnuplot.Register(interfaced);
OL.client Gnuplot/InputFiles.Set(f.plt);
OL.iftrue(TRANSIENT)
	OL.client PostArray.List(fmax.txt,1,2,Solution/tmin, fmax.txt,1,8,Solution/fmin);
OL.else
	OL.client Gnuplot.Active(0); // no gnuplot
	OL.client PostArray.List(fmax.txt,1,8,Solution/fobj);
OL.endif

