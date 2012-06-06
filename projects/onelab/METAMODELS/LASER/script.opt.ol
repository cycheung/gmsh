Plugin(MathEval).Expression0= "v0";
Plugin(MathEval).Expression1= "";
Plugin(MathEval).Expression2= "";
Plugin(MathEval).Expression3= "";
Plugin(MathEval).Expression4= "";
Plugin(MathEval).Expression5= "";
Plugin(MathEval).Expression6= "";
Plugin(MathEval).Expression7= "";
Plugin(MathEval).Expression8= "";
Plugin(MathEval).TimeStep=OL.eval(OL.get(Parameters/Elmer/NumStep)-1);
Plugin(MathEval).View=0;
Plugin(MathEval).OtherTimeStep=-1;
Plugin(MathEval).OtherView=-1;
Plugin(MathEval).ForceInterpolation=0;
Plugin(MathEval).PhysicalRegion=-1;
Plugin(MathEval).Run;

Plugin(CutPlane).A=0;
Plugin(CutPlane).B=-1;
Plugin(CutPlane).C=0;
Plugin(CutPlane).D= OL.get(PostPro/ZSURF0);
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=1;
Plugin(CutPlane).Run;

Save View [2] "tempsurf.txt";

Plugin(MinMax).View=2;
Plugin(MinMax).OverTime=0;
Plugin(MinMax).Argument=0;
Plugin(MinMax).Run;

Save View [3] "tempmin.txt";
Save View [4] "tempmax.txt";
