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
Plugin(CutPlane).D= OL.get(PostPro/ZSURF0);
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

//CUT THE Z_PLANES
For k In {1:5}
    Plugin(CutPlane).A=0;
    Plugin(CutPlane).B=-1;
    Plugin(CutPlane).C=0;
If(k==1)
    Plugin(CutPlane).D= OL.get(PostPro/ZSURF0);
EndIf
If(k==2)
    Plugin(CutPlane).D= OL.get(PostPro/ZSURF1);
EndIf
If(k==3)
    Plugin(CutPlane).D= OL.get(PostPro/ZSURF2);
EndIf
If(k==4)
    Plugin(CutPlane).D= OL.get(PostPro/ZSURF3);
EndIf
If(k==5)
    Plugin(CutPlane).D= OL.get(PostPro/ZSURF4);
EndIf
    Plugin(CutPlane).ExtractVolume=0;
    Plugin(CutPlane).RecurLevel=4;
    Plugin(CutPlane).TargetError=0;
    Plugin(CutPlane).View=0;
    Plugin(CutPlane).Run;
    Save View [4+k] Sprintf("temp%g.txt", k-1);
EndFor

//MATHEVAL THE CUTPLANE WITH STEP
For k In {1:5}
 Plugin(MathEval).Expression0= "Step(v0-320.15)*2*pi*x";
 Plugin(MathEval).TimeStep=-1;
 Plugin(MathEval).View=4+k;
 Plugin(MathEval).OtherTimeStep=-1;
 Plugin(MathEval).OtherView=-1;
 Plugin(MathEval).ForceInterpolation=0;
 Plugin(MathEval).PhysicalRegion=-1;
 Plugin(MathEval).Run;
 Save View [9+k] Sprintf("step%g.txt", k-1);
EndFor

//INTEGRATE
For k In {1:5}
 Plugin(Integrate).View=9+k;
 Plugin(Integrate).Run; 
 Save View [14+k] Sprintf("active%g.txt", k-1);
EndFor