ListDepth=OL.get(PostPro/4ZSURF,choices.expand( '{comma}' ));
nbDepth = #ListDepth[];
Printf("NbDepth = %g", nbDepth);

NB = 21;
SW = OL.get(PostPro/3SKINWIDTH)-1e-4;
DZ = 0.2 / NB;

For k In {0:NB}
  ListDepthF[k] = SW - k*DZ;
EndFor

nbDepthF = #ListDepthF[];
Printf("The %g considered depths are:", nbDepthF);
For k In {0:NB}
  Printf("%f", ListDepthF[k]);
EndFor

// EXTRACT TEMPERATURE FIELD AT PROBETIME
Plugin(ExtractElements).MinVal=0;
Plugin(ExtractElements).MaxVal=0;
Plugin(ExtractElements).TimeStep = OL.eval(abs(OL.get(PostPro/2PROBETIME)/(OL.get(Parameters/Elmer/2TimeStep)*1e-3)));
Plugin(ExtractElements).Visible=1;
Plugin(ExtractElements).Dimension=2;
Plugin(ExtractElements).View=0;
Plugin(ExtractElements).Run;

//EVALUATE TEMPERATURE PROFILE AT: SURFACE, INTERFACE AND PROBE DEPTH
// View 1 -> nbDepth
For k In {1:nbDepth}
 Plugin(CutPlane).A=0;
 Plugin(CutPlane).B=-1;
 Plugin(CutPlane).C=0;
 Plugin(CutPlane).D= ListDepth[k-1]*1e-3;
 Plugin(CutPlane).ExtractVolume=0;
 Plugin(CutPlane).RecurLevel=4;
 Plugin(CutPlane).TargetError=0;
 Plugin(CutPlane).View=1;
 Plugin(CutPlane).Run;
 Save View [1+k] Sprintf("templaser%g.txt", k-1);
EndFor

ViewNum=nbDepth+1;

//EVALUATE TEMPERATURE PROFILE AT NB EQUIDISTANT DEPTHS
For k In {1:nbDepthF}
    Plugin(CutPlane).A=0;
    Plugin(CutPlane).B=-1;
    Plugin(CutPlane).C=0;
    Plugin(CutPlane).D= ListDepthF[k-1]*1e-3;
    Plugin(CutPlane).ExtractVolume=0;
    Plugin(CutPlane).RecurLevel=4;
    Plugin(CutPlane).TargetError=0;
    Plugin(CutPlane).View=1;
    Plugin(CutPlane).Run;
    Save View [ViewNum+k] Sprintf("temp%g.txt", k-1);
EndFor

//EVALUATE FOR EACH DEPTH THE BOOLEAN FUNCTION: ABOVE THRESHOLD (1/0)
//MULTIPLY BY 2 PI X IN VIEW OF INTEGRATION
For k In {1:nbDepthF}
 Plugin(MathEval).Expression0= "Step(v0-OL.get(PostPro/5OVERTEMP))*2*pi*x";
 Plugin(MathEval).TimeStep=-1;
 Plugin(MathEval).View=ViewNum+k;
 Plugin(MathEval).OtherTimeStep=-1;
 Plugin(MathEval).OtherView=-1;
 Plugin(MathEval).ForceInterpolation=0;
 Plugin(MathEval).PhysicalRegion=-1;
 Plugin(MathEval).Run;
EndFor

ViewNum=ViewNum+nbDepthF;

//INTEGRATE
For k In {1:nbDepthF}
 Plugin(Integrate).View=ViewNum+k;
 Plugin(Integrate).Run; 
EndFor

ViewNum=ViewNum+nbDepthF;

//MIN MAX 
For k In {1:nbDepthF}
 Plugin(MinMax).View=ViewNum+k;
 Plugin(MinMax).OverTime=1;
 Plugin(MinMax).Argument=0;
 Plugin(MinMax).Run;
 Save View [ViewNum+nbDepthF+(k*2)] Sprintf("activeMax%g.txt", k-1);
EndFor

ViewNum=ViewNum+nbDepthF;

Combine ElementsByViewName;
Save View [6] "activeMax.txt"; 

