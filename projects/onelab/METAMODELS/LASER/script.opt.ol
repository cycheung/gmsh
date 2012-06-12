//SURF TEMP AT THE END OF LASER APPLICATION
Plugin(MathEval).Expression0= "v0";
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
Plugin(CutPlane).D= OL.get(PostPro/ZSURF,choices.comp(0));
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=1;
Plugin(CutPlane).Run;
Save View [2] "tempsurf.txt";

//MINMAX SURF TEMP
Plugin(MinMax).View=2;
Plugin(MinMax).OverTime=0;
Plugin(MinMax).Argument=0;
Plugin(MinMax).Run;
Save View [3] "tempmin.txt";
Save View [4] "tempmax.txt";

ListDepth=OL.get(PostPro/ZSURF,choices.expand( '{comma}' ));

ViewNum=4;

//CUT THE Z_PLANES
For k In {1:# ListDepth[]}
    Plugin(CutPlane).A=0;
    Plugin(CutPlane).B=-1;
    Plugin(CutPlane).C=0;
    Plugin(CutPlane).D= ListDepth[k];
    Plugin(CutPlane).ExtractVolume=0;
    Plugin(CutPlane).RecurLevel=4;
    Plugin(CutPlane).TargetError=0;
    Plugin(CutPlane).View=0;
    Plugin(CutPlane).Run;
    Save View [ViewNum+k] Sprintf("temp%g.txt", k-1);
EndFor

ViewNum=ViewNum+# ListDepth[];

//MATHEVAL THE CUTPLANE WITH STEP
For k In {1:# ListDepth[]}
 Plugin(MathEval).Expression0= "Step(v0-320.15)*2*pi*x";
 Plugin(MathEval).TimeStep=-1;
 Plugin(MathEval).View=4+k;
 Plugin(MathEval).OtherTimeStep=-1;
 Plugin(MathEval).OtherView=-1;
 Plugin(MathEval).ForceInterpolation=0;
 Plugin(MathEval).PhysicalRegion=-1;
 Plugin(MathEval).Run;
 //Save View [ViewNum+k] Sprintf("step%g.txt", k-1);
EndFor

ViewNum=ViewNum+# ListDepth[];

//INTEGRATE
For k In {1:# ListDepth[]}
 Plugin(Integrate).View=9+k;
 Plugin(Integrate).Run; 
 Save View [ViewNum+k] Sprintf("active%g.txt", k-1);
EndFor

ViewNum=ViewNum+# ListDepth[];

//MIN MAX 
For k In {1:# ListDepth[]}
 Plugin(MinMax).View=14+k;
 Plugin(MinMax).OverTime=1;
 Plugin(MinMax).Argument=0;
 Plugin(MinMax).Run;
 Save View [ViewNum+(k*2)] Sprintf("activeMax%g.txt", k-1);
EndFor

ViewNum=ViewNum+# ListDepth[];

Combine ElementsByViewName;
Save View [9] "activeMax.txt";  // toujours View[9] ??