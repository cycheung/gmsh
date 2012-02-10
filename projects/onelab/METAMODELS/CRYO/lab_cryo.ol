OL.client Mesher.Register(encapsulated);
OL.client Mesher/InputFiles.Set(OL.getValue(Arguments/FileName).geo);
%OL.client Mesher/LineOptions.Set(-v 0);

OL.client ElmerGrid.Register(interfaced);
OL.client ElmerGrid/PreLineOptions.Set(14 2);
OL.client ElmerGrid/InputFiles.Set(OL.getValue(Arguments/FileName).msh);
OL.client ElmerGrid/LineOptions.Set(-out mesh);
OL.client ElmerGrid/OutputFiles.Set(mesh/mesh.header);

OL.client Elmer.Register(interfaced);
OL.client Elmer/InputFiles.Set(OL.getValue(Arguments/FileName).sif.ol);

OL.client Post.Register(interfaced);
OL.client Post/InputFiles.Set(solution.pos,script.opt);
OL.client Post/LineOptions.Set("-");

OL.client Gnuplot.Register(interfaced);
OL.client Gnuplot/InputFiles.Set(plot.gnu);

OL.client Matlab.Register(interfaced);
OL.client Matlab/InputFiles.Set(plotMatlab.m);
OL.client Matlab/LineOptions.Set(-nosplash -r plotMatlab);

OL.client PostArray.List(fmax.txt,1,2,Solution/tmin, fmax.txt,1,8,Solution/fmin);

OL.client GmshMerge/InputFiles.Set(solution.pos);