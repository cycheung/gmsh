ListDepth=OL.get(PostPro/ZSURF,choices.expand( '{comma}' ));
nbDepth = #ListDepth[];

OL.block
 SKINWIDTH.number(, Parameters/Skin/,"Skin width"); 
 SKINWIDTH.setValue(OL.eval((OL.get(Parameters/Skin/3DERMIS)+OL.get(Parameters/Skin/2EPIDERMIS))*1e-3));
 SKINWIDTH.setVisible(0);
 ZSURFF.number( , PostPro/,"Z full coordinates"); 
 ZSURFF.resetChoices();
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.00099 * 1e-3) );  
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.00601 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.01249 * 1e-3) );  
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.01801 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.02599 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.03749 * 1e-3) );  
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.04349 * 1e-3) ); 	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.05009 * 1e-3) );  
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.05609 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.06249 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.07501 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.08751 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.10001 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.11251 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.12501 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.13751 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.15001 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.16251 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.17501 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.18751 * 1e-3) );  	 	 
 ZSURFF.addChoices( OL.eval( OL.get(SKINWIDTH) - 0.20001 * 1e-3) );
 ZSURFF.setVisible(0);
OL.endblock

ListDepthF=OL.get(PostPro/ZSURFF,choices.expand( '{comma}' ));
nbDepthF = #ListDepthF[];

Plugin(ExtractElements).MinVal=0;
Plugin(ExtractElements).MaxVal=0;
Plugin(ExtractElements).TimeStep = OL.eval(OL.get(PostPro/PROBETIME)/OL.get(Parameters/Elmer/TimeStep)*1000-1);
Plugin(ExtractElements).Visible=1;
Plugin(ExtractElements).Dimension=2;
Plugin(ExtractElements).View=0;
Plugin(ExtractElements).Run;

//SKINTEMP AT THE END OF LASER APPLICATION
// View 1 -> nbDepth
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

ViewNum=nbDepth+1;

//CUT THE Z_PLANES
For k In {1:nbDepthF}
    Plugin(CutPlane).A=0;
    Plugin(CutPlane).B=-1;
    Plugin(CutPlane).C=0;
    Plugin(CutPlane).D= ListDepthF[k-1];
    Plugin(CutPlane).ExtractVolume=0;
    Plugin(CutPlane).RecurLevel=4;
    Plugin(CutPlane).TargetError=0;
    Plugin(CutPlane).View=1;
    Plugin(CutPlane).Run;
    Save View [ViewNum+k] Sprintf("temp%g.txt", k-1);
EndFor

//MATHEVAL THE CUTPLANE WITH STEP
For k In {1:nbDepthF}
 Plugin(MathEval).Expression0= "Step(v0-OL.get(Parameters/Skin/OVERTEMP)+1)*2*pi*x";
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

