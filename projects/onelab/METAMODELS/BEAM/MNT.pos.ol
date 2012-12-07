OL.block
CUTPLANE.radioButton(1,,"Use CutPlane instead of CutGrid");
CUTPLANE.setVisible(0);
OL.endblock

k=PostProcessing.NbViews; 
Printf("%g views in MNT.pos (should be 9)", k);

// Calculation of M
// View[k]

OL.iftrue(CUTPLANE)
Plugin(CutPlane).A=-1;
Plugin(CutPlane).B=0;
Plugin(CutPlane).C=0;
Plugin(CutPlane).D= OL.get(1Geometry/X);
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=k-9;
Plugin(CutPlane).Run;
OL.else
Plugin(CutGrid).X0=OL.get(1Geometry/X);
Plugin(CutGrid).Y0=-OL.eval(0.5 * OL.get(1Geometry/A) );
Plugin(CutGrid).Z0=-OL.eval(0.5 * OL.get(1Geometry/B) );
Plugin(CutGrid).X1=OL.get(1Geometry/X);
Plugin(CutGrid).Y1=-OL.eval(0.5 * OL.get(1Geometry/A));
Plugin(CutGrid).Z1=OL.eval(0.5 * OL.get(1Geometry/B));
Plugin(CutGrid).X2=OL.get(1Geometry/X);
Plugin(CutGrid).Y2=OL.eval(0.5 * OL.get(1Geometry/A));
Plugin(CutGrid).Z2=-OL.eval(0.5 * OL.get(1Geometry/B));
Plugin(CutGrid).NumPointsU=100;
Plugin(CutGrid).NumPointsV=100;
Plugin(CutGrid).ConnectPoints=1;
Plugin(CutGrid).View=k-9;
Plugin(CutGrid).Run;
OL.endif

//View[k+1]
Plugin(MathEval).Expression0= "v0*y";
Plugin(MathEval).Expression1= "";
Plugin(MathEval).Expression2= "";
Plugin(MathEval).Expression3= "";
Plugin(MathEval).Expression4= "";
Plugin(MathEval).Expression5= "";
Plugin(MathEval).Expression6= "";
Plugin(MathEval).Expression7= "";
Plugin(MathEval).Expression8= "";
Plugin(MathEval).TimeStep=-1;
Plugin(MathEval).View=k;
Plugin(MathEval).OtherTimeStep=-1;
Plugin(MathEval).OtherView=-1;
Plugin(MathEval).ForceInterpolation=0;
Plugin(MathEval).PhysicalRegion=-1;
Plugin(MathEval).Run;

//View[k+2]
Plugin(Integrate).View=k+1;
Plugin(Integrate).OverTime=-1;
Plugin(Integrate).Dimension=2;
Plugin(Integrate).Run;

// Calculation of N
//View[k+3]

OL.iftrue(CUTPLANE)
Plugin(CutPlane).A=-1;
Plugin(CutPlane).B=0;
Plugin(CutPlane).C=0;
Plugin(CutPlane).D= OL.get(1Geometry/X);
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=k-9;
Plugin(CutPlane).Run;
OL.else
Plugin(CutGrid).X0=OL.get(1Geometry/X);
Plugin(CutGrid).Y0=-OL.eval(0.5 * OL.get(1Geometry/A));
Plugin(CutGrid).Z0=-OL.eval(0.5 * OL.get(1Geometry/B));
Plugin(CutGrid).X1=OL.get(1Geometry/X);
Plugin(CutGrid).Y1=-OL.eval(0.5 * OL.get(1Geometry/A));
Plugin(CutGrid).Z1=OL.eval(0.5 * OL.get(1Geometry/B));
Plugin(CutGrid).X2=OL.get(1Geometry/X);
Plugin(CutGrid).Y2=OL.eval(0.5 * OL.get(1Geometry/A));
Plugin(CutGrid).Z2=-OL.eval(0.5 * OL.get(1Geometry/B));
Plugin(CutGrid).NumPointsU=100;
Plugin(CutGrid).NumPointsV=100;
Plugin(CutGrid).ConnectPoints=1;
Plugin(CutGrid).View=k-9;
Plugin(CutGrid).Run;
OL.endif


//View[k+4]
Plugin(Integrate).View=k+3;
Plugin(Integrate).OverTime=-1;
Plugin(Integrate).Dimension=2;
Plugin(Integrate).Run;

// Calculation of T
//View[k+5]
OL.iftrue(CUTPLANE)
Plugin(CutPlane).A=-1;
Plugin(CutPlane).B=0;
Plugin(CutPlane).C=0;
Plugin(CutPlane).D= OL.get(1Geometry/X);
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=k-6;
Plugin(CutPlane).Run;
OL.else
Plugin(CutGrid).X0=OL.get(1Geometry/X);
Plugin(CutGrid).Y0=-OL.eval(0.5 * OL.get(1Geometry/A));
Plugin(CutGrid).Z0=-OL.eval(0.5 * OL.get(1Geometry/B));
Plugin(CutGrid).X1=OL.get(1Geometry/X);
Plugin(CutGrid).Y1=-OL.eval(0.5 * OL.get(1Geometry/A));
Plugin(CutGrid).Z1=OL.eval(0.5 * OL.get(1Geometry/B));
Plugin(CutGrid).X2=OL.get(1Geometry/X);
Plugin(CutGrid).Y2=OL.eval(0.5 * OL.get(1Geometry/A));
Plugin(CutGrid).Z2=-OL.eval(0.5 * OL.get(1Geometry/B));
Plugin(CutGrid).NumPointsU=100;
Plugin(CutGrid).NumPointsV=100;
Plugin(CutGrid).ConnectPoints=1;
Plugin(CutGrid).View=k-6;
Plugin(CutGrid).Run;
OL.endif


//View[k+6]
Plugin(Integrate).View=k+5;
Plugin(Integrate).OverTime=-1;
Plugin(Integrate).Dimension=2;
Plugin(Integrate).Run;

//View[k+7]
Plugin(Scal2Vec).NameNewView= "MNT";
Plugin(Scal2Vec).ViewX=k+2;
Plugin(Scal2Vec).ViewY=k+4;
Plugin(Scal2Vec).ViewZ=k+6;
Plugin(Scal2Vec).Run;

Save View [k+7] "MNT.txt"; 

