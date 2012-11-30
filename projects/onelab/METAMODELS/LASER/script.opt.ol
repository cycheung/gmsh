ListDepth=OL.get(PostPro/ZSURF,choices.expand( '{comma}' ));
nbDepth = #ListDepth[];

//SKINTEMP AT THE END OF LASER APPLICATION
Plugin(MathEval).Expression0= "v0";
Plugin(MathEval).TimeStep=OL.eval(OL.get(Parameters/Elmer/NumStep)-1);
Plugin(MathEval).View=0;
Plugin(MathEval).OtherTimeStep=-1;
Plugin(MathEval).OtherView=-1;
Plugin(MathEval).ForceInterpolation=0;
Plugin(MathEval).PhysicalRegion=-1;
Plugin(MathEval).Run;

For k In {1:nbDepth}
 Plugin(CutPlane).A=0;
 Plugin(CutPlane).B=-1;
 Plugin(CutPlane).C=0;
 Plugin(CutPlane).D= ListDepth[k-1];
 Plugin(CutPlane).ExtractVolume=0;
 Plugin(CutPlane).RecurLevel=4;
 Plugin(CutPlane).TargetError=0;
 Plugin(CutPlane).View=1;
 Plugin(CutPlane).Run;
 Save View [1+k] Sprintf("templaser%g.txt", k-1);
EndFor

//MINMAX SURF TEMP
Plugin(MinMax).View=2;
Plugin(MinMax).OverTime=0;
Plugin(MinMax).Argument=0;
Plugin(MinMax).Run;
//Save View [nbDepth+2] "tempmin.txt";
//Save View [nbDepth+3] "tempmax.txt";

ViewNum=nbDepth+3;

//CUT THE Z_PLANES
For k In {1:nbDepth}
    Plugin(CutPlane).A=0;
    Plugin(CutPlane).B=-1;
    Plugin(CutPlane).C=0;
    Plugin(CutPlane).D= ListDepth[k-1];
    Plugin(CutPlane).ExtractVolume=0;
    Plugin(CutPlane).RecurLevel=4;
    Plugin(CutPlane).TargetError=0;
    Plugin(CutPlane).View=0;
    Plugin(CutPlane).Run;
    Save View [ViewNum+k] Sprintf("temp%g.txt", k-1);
EndFor

//MATHEVAL THE CUTPLANE WITH STEP
For k In {1:nbDepth}
 Plugin(MathEval).Expression0= "Step(v0-320.15)*2*pi*x";
 Plugin(MathEval).TimeStep=-1;
 Plugin(MathEval).View=ViewNum+k;
 Plugin(MathEval).OtherTimeStep=-1;
 Plugin(MathEval).OtherView=-1;
 Plugin(MathEval).ForceInterpolation=0;
 Plugin(MathEval).PhysicalRegion=-1;
 Plugin(MathEval).Run;
EndFor

ViewNum=ViewNum+nbDepth;

//INTEGRATE
For k In {1:nbDepth}
 Plugin(Integrate).View=ViewNum+k;
 Plugin(Integrate).Run; 
EndFor

ViewNum=ViewNum+nbDepth;

//MIN MAX 
For k In {1:nbDepth}
 Plugin(MinMax).View=ViewNum+k;
 Plugin(MinMax).OverTime=1;
 Plugin(MinMax).Argument=0;
 Plugin(MinMax).Run;
 Save View [ViewNum+nbDepth+(k*2)] Sprintf("activeMax%g.txt", k-1);
EndFor

ViewNum=ViewNum+nbDepth;

Combine ElementsByViewName;
Save View [9] "activeMax.txt";  // toujours View[9] ??