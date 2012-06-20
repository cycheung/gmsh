Plugin(MathEval).Expression0= "Step(v0-OL.get(Parameters/Model/OVERTEMP))";
Plugin(MathEval).TimeStep=-1;
Plugin(MathEval).View=0;
Plugin(MathEval).OtherTimeStep=-1;
Plugin(MathEval).OtherView=-1;
Plugin(MathEval).ForceInterpolation=0;
Plugin(MathEval).PhysicalRegion=-1;
Plugin(MathEval).Run;
Save View [1] "overheat.pos";

Plugin(CutPlane).A=-1;
Plugin(CutPlane).B=0;
Plugin(CutPlane).C=0;
Plugin(CutPlane).D=0;
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=1;
Plugin(CutPlane).Run;
Save View [2] "step.txt";

Plugin(Integrate).View=2;
Plugin(Integrate).OverTime=1;
Plugin(Integrate).Run; 
Save View [3] "duration.txt";
