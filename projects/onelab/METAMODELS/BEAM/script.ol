k=PostProcessing.NbViews; 

// Calculation of M

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
Plugin(MathEval).View=k-12;
Plugin(MathEval).OtherTimeStep=-1;
Plugin(MathEval).OtherView=-1;
Plugin(MathEval).ForceInterpolation=0;
Plugin(MathEval).PhysicalRegion=-1;
Plugin(MathEval).Run;

Plugin(CutPlane).A=-1;
Plugin(CutPlane).B=0;
Plugin(CutPlane).C=0;
Plugin(CutPlane).D= OL.get(1Geometry/X);
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=k;
Plugin(CutPlane).Run;

Plugin(Integrate).View=k+1;
Plugin(Integrate).OverTime=-1;
Plugin(Integrate).Run;

// Calculation of N

Plugin(CutPlane).A=-1;
Plugin(CutPlane).B=0;
Plugin(CutPlane).C=0;
Plugin(CutPlane).D= OL.get(1Geometry/X);
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=k-9;
Plugin(CutPlane).Run;

Plugin(Integrate).View=k+3;
Plugin(Integrate).OverTime=-1;
Plugin(Integrate).Run;

// Calculation of T

Plugin(CutPlane).A=-1;
Plugin(CutPlane).B=0;
Plugin(CutPlane).C=0;
Plugin(CutPlane).D= OL.get(1Geometry/X);
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=k-6;
Plugin(CutPlane).Run;

Plugin(Integrate).View=k+5;
Plugin(Integrate).OverTime=-1;
Plugin(Integrate).Run;

// 

Plugin(Scal2Vec).NameNewView= "MNT";
Plugin(Scal2Vec).ViewX=k+2;
Plugin(Scal2Vec).ViewY=k+4;
Plugin(Scal2Vec).ViewZ=k+6;
Plugin(Scal2Vec).Run;

Save View [k+7] "MNT.txt"; 

