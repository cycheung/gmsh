% global variables 

OL.parameter TRANSIENT.number(0,,Flag transient simulation); TRANSIENT.AddChoices(0,1);
OL.parameter DISTANT.number(0,,Flag remote host); DISTANT.AddChoices(0,1);

% DISTANT=1 : haddock.mema.ucl.ac.be  130.104.237.29
%DISTANT=2: localhost

OL.ifequal(DISTANT,1)
    OL.parameter HOST.string( francois@haddock.mema.ucl.ac.be ); 
    OL.parameter WORKDIR.string( WORK/FH/CRYO );
    OL.parameter BINDIR.string( /home/francois/bin );
OL.endif

OL.ifequal(DISTANT,2)
% The remote client can be localhost.
% This allows solving on the local machine but in another directory.
    OL.parameter HOST.string( frhenrotte@localhost ); 
    OL.parameter WORKDIR.string( tmp/CRYO );
    OL.parameter BINDIR.string( ~/SCRIPTS/bash );
OL.endif

%1) Client Gmsh pour le maillage initial
OL.iftrue(DISTANT)
    OL.client Mesher.Register(encapsulated,OL.getValue(HOST),OL.getValue(WORKDIR),OL.getValue(BINDIR)/gmsh);
OL.else
    OL.client Mesher.Register(encapsulated,gmsh);
OL.endif
OL.client Mesher.In( ./OL.getValue(Arguments/FileName).geo);
OL.client Mesher.Args(OL.getValue(Arguments/FileName).geo);
OL.client Mesher.Out(OL.getValue(Arguments/FileName).msh);

%2) Client ElmerGrid pour convertir le maillage pour Elmer
OL.iftrue(DISTANT)
    OL.client ElmerGrid.Register(interfaced,OL.getValue(HOST),OL.getValue(WORKDIR),OL.getValue(BINDIR)/ElmerGrid);
OL.else
    OL.client ElmerGrid.Register(interfaced,ElmerGrid);
OL.endif
OL.client ElmerGrid.In( OL.getValue(Arguments/FileName).msh);
OL.client ElmerGrid.Args(14 2 OL.getValue(Arguments/FileName).msh -out mesh);
OL.client ElmerGrid.Out(mesh);

%3) Client Elmer pour lancer la simulation avec Elmer
OL.iftrue(DISTANT)
    OL.client Elmer.Register(interfaced,OL.getValue(HOST),OL.getValue(WORKDIR),OL.getValue(BINDIR)/ElmerSolver);
OL.else
    OL.client Elmer.Register(interfaced,ElmerSolver);
OL.endif
OL.client Elmer.In( ./ELMERSOLVER_STARTINFO.ol , ./OL.getValue(Arguments/FileName).sif.ol);
OL.client Elmer.Out( solution.pos , tempevol.txt);

%4)Client Postpro pour executer un script gmsh
OL.iftrue(DISTANT)
    OL.client Post.Register(interfaced,OL.getValue(HOST),OL.getValue(WORKDIR),OL.getValue(BINDIR)/gmsh);
OL.else
    OL.client Post.Register(interfaced,gmsh);
OL.endif
OL.client Post.In( solution.pos , ./script_cryo.opt );
OL.client Post.Args(solution.pos script_cryo.opt -);
OL.client Post.Out( ./f.txt , ./fmax.txt);

%5)Client Postpro pour charger la solution calcul avec Elmer 
OL.ifntrue(DISTANT)
    OL.client GmshMerge.In(solution.pos);
OL.endif

%6)Client Postpro pour extraire la min et la max de la courbe
OL.client Gnuplot.Register(interfaced,gnuplot);
OL.client Gnuplot.In(f.plt);
OL.iftrue(TRANSIENT)
        OL.client PostArray.List(fmax.txt,-1,2,Solution/tmin, fmax.txt,-1,8,Solution/fmin);
OL.else
        OL.client Gnuplot.Active(0); // no gnuplot
        OL.client PostArray.List(aaa,-1,8,Solution/fobj);
OL.endif
