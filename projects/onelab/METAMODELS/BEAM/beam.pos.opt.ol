k=PostProcessing.NbViews; 

//Printf("%g views in beam.pos.opt", k);
//For ii In {10:k}
//Delete View[k-ii]; 
//EndFor

View[k-9].Visible = 0;
View[k-8].Visible = 0;
View[k-7].Visible = 0;
View[k-6].Visible = 0;
View[k-5].Visible = 0;
View[k-4].Visible = 0;
View[k-3].Visible = 0;
View[k-2].ShowElement = 1;
View[k-2].VectorType = 5;
View[k-2].ExternalView = k-9;
OL.block
RANGE.number();
RANGE.setValue(OL.eval(fabs((OL.get(9Results/vmax))-(OL.get(9Results/vmin))) ));
RANGE.setVisible(0);
OL.endblock
OL.if( OL.get(RANGE) < 1e-7 ) 
  View[k-2].DisplacementFactor = 10000 ;
OL.else
  View[k-2].DisplacementFactor = OL.eval(0.33*OL.get(1Geometry/L)/OL.get(RANGE) ) ;
OL.endif
View[k-1].Visible = 0;

Plugin(CutGrid).X0=OL.get(1Geometry/X);
Plugin(CutGrid).Y0=OL.eval(-OL.get(1Geometry/A)*0.5);
Plugin(CutGrid).Z0=OL.eval(-OL.get(1Geometry/B)*0.5);
Plugin(CutGrid).X1=OL.get(1Geometry/X);
Plugin(CutGrid).Y1=OL.eval( OL.get(1Geometry/A)*0.5);
Plugin(CutGrid).Z1=OL.eval(-OL.get(1Geometry/B)*0.5);
Plugin(CutGrid).X2=OL.get(1Geometry/X);
Plugin(CutGrid).Y2=OL.eval(-OL.get(1Geometry/A)*0.5);
Plugin(CutGrid).Z2=OL.eval( OL.get(1Geometry/B)*0.5);
Plugin(CutGrid).NumPointsU=10;
Plugin(CutGrid).NumPointsV=10;
Plugin(CutGrid).ConnectPoints=1;
Plugin(CutGrid).View=k-9;
Plugin(CutGrid).Run;

Plugin(CutGrid).X0=OL.get(1Geometry/X);
Plugin(CutGrid).Y0=OL.eval(-OL.get(1Geometry/A)*0.5);
Plugin(CutGrid).Z0=OL.eval(-OL.get(1Geometry/B)*0.5);
Plugin(CutGrid).X1=OL.get(1Geometry/X);
Plugin(CutGrid).Y1=OL.eval( OL.get(1Geometry/A)*0.5);
Plugin(CutGrid).Z1=OL.eval(-OL.get(1Geometry/B)*0.5);
Plugin(CutGrid).X2=OL.get(1Geometry/X);
Plugin(CutGrid).Y2=OL.eval(-OL.get(1Geometry/A)*0.5);
Plugin(CutGrid).Z2=OL.eval( OL.get(1Geometry/B)*0.5);
Plugin(CutGrid).NumPointsU=10;
Plugin(CutGrid).NumPointsV=10;
Plugin(CutGrid).ConnectPoints=1;
Plugin(CutGrid).View=k-6;
Plugin(CutGrid).Run;

View[k].Visible = 0;
View[k+1].Visible = 0;

Plugin(Scal2Vec).NameNewView= "Tx X=OL.get(1Geometry/X)";
Plugin(Scal2Vec).ViewX=k;
Plugin(Scal2Vec).ViewY=k+1;
Plugin(Scal2Vec).ViewZ=-1;
Plugin(Scal2Vec).Run;

View[k+2].VectorType = 3;
View[k+2].ShowScale = 0;
View[k+2].Visible = 1;