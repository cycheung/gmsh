// View 1
Plugin(MathEval).Expression0= "Step(v0-OL.get(Parameters/Skin/OVERTEMP))";
Plugin(MathEval).TimeStep=-1;
Plugin(MathEval).View=0;
Plugin(MathEval).OtherTimeStep=-1;
Plugin(MathEval).OtherView=-1;
Plugin(MathEval).ForceInterpolation=0;
Plugin(MathEval).PhysicalRegion=-1;
Plugin(MathEval).Run;
Save View [1] "aboveThres.pos";

// View 2
Plugin(CutPlane).A=-1;
Plugin(CutPlane).B=0;
Plugin(CutPlane).C=0;
Plugin(CutPlane).D=0;
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=1;
Plugin(CutPlane).Run;
//Save View [2] "step.txt";

// View 3
Plugin(Integrate).View=2;
Plugin(Integrate).OverTime=1;
Plugin(Integrate).Dimension=1;
Plugin(Integrate).Run;
Save View [3] "duration.txt";

// View 4
Plugin(MathEval).Expression0= "Step(v0-OL.get(Parameters/Skin/OVERTEMP))*6.2831853e9*x";
Plugin(MathEval).TimeStep=-1;
Plugin(MathEval).View=0;
Plugin(MathEval).OtherTimeStep=-1;
Plugin(MathEval).OtherView=-1;
Plugin(MathEval).ForceInterpolation=0;
Plugin(MathEval).PhysicalRegion=-1;
Plugin(MathEval).Run;

// View 5
Plugin(Integrate).View=4;
Plugin(Integrate).OverTime=-1;
Plugin(Integrate).Dimension=2;
Plugin(Integrate).Run; 

Save View [5] "volume.txt" ;

// View 6 and 7
Plugin(MinMax).View=5;
Plugin(MinMax).OverTime=1;
Plugin(MinMax).Argument=0;
Plugin(MinMax).Run;
Save View [7] "volmax.txt";

// View 8
Plugin(MathEval).Expression0= "Step(v0-OL.get(Parameters/Skin/OVERTEMP))*6.2831853e9*x";
Plugin(MathEval).TimeStep=-1;
Plugin(MathEval).View=0;
Plugin(MathEval).OtherTimeStep=-1;
Plugin(MathEval).OtherView=-1;
Plugin(MathEval).ForceInterpolation=0;
Plugin(MathEval).PhysicalRegion=1;
Plugin(MathEval).Run;

// View 9
Plugin(Integrate).View=8;
Plugin(Integrate).OverTime=-1;
Plugin(Integrate).Dimension=2;
Plugin(Integrate).Run; 

Save View [9] "volderm.txt" ;

// View 10 and 11
Plugin(MinMax).View=9;
Plugin(MinMax).OverTime=1;
Plugin(MinMax).Argument=0;
Plugin(MinMax).Run;
Save View [11] "voldermmax.txt";