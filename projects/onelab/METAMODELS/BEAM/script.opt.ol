

k=PostProcessing.NbViews; 
Printf("%g views in script.pos", k);

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

#OL.block
#RANGE.number();
#RANGE.setValue(OL.eval(fabs((OL.get(9Results/vmax))-(OL.get(9Results/vmin))) ));
#RANGE.setVisible(0);
#OL.endblock
#OL.if( OL.get(RANGE) < 1e-7 ) 
#  View[k-2].DisplacementFactor = 10000 ;
#OL.else
#  View[k-2].DisplacementFactor = OL.eval(0.33*OL.get(1Geometry/L)/OL.get(RANGE) ) ;
#OL.endif

OL.block
MAGN.number(2000, 1Geometry/,"Magnification factor");
OL.endblock

View[k-2].DisplacementFactor = OL.get(MAGN) ;

View[k-1].Visible = 0;

//View[k]
Plugin(CutPlane).A=-1;
Plugin(CutPlane).B=0;
Plugin(CutPlane).C=0;
Plugin(CutPlane).D= OL.get(1Geometry/X);
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=k-9;
Plugin(CutPlane).Run;

//View[k+1]
Plugin(CutPlane).A=-1;
Plugin(CutPlane).B=0;
Plugin(CutPlane).C=0;
Plugin(CutPlane).D= OL.get(1Geometry/X);
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=k-6;
Plugin(CutPlane).Run;

View[k].Visible = 0;
View[k+1].Visible = 0;

Printf("sigma_xx Min = %g Max = %g", View[k].Min, View[k].Max);
Printf("sigma_xy Min = %g Max = %g", View[k+1].Min, View[k+1].Max);

//View[k+2]
Plugin(Scal2Vec).NameNewView= "Tx X=OL.get(1Geometry/X)";
Plugin(Scal2Vec).ViewX=k;
Plugin(Scal2Vec).ViewY=k+1;
Plugin(Scal2Vec).ViewZ=-1;
Plugin(Scal2Vec).Run;

OL.iftrue(LOOP)
View[k+2].RangeType = 2;
View[k+2].CustomMin = -1e6;
View[k+2].CustomMax = 1e6;
View[k+2].ShowScale = 0;
OL.endif
View[k+2].Visible = 1;
View[k+2].VectorType = 3;

