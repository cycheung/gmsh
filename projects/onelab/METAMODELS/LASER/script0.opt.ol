
// computation of the scalar output values of the metamodel

// MAX TEMP AT DERMIS-EPIDERMIS INTERFACE

// DefineConstant [DERMIS = {0.05, Name "Parameters/Skin/3Mesh size"}];

Plugin(Probe).View=0;
Plugin(Probe).X=1.e-6;
Plugin(Probe).Y=(OL.get(Parameters/Skin/3DERMIS)-1e-4)*1.e-3;
Plugin(Probe).Z=0;
Plugin(Probe).Run;

Plugin(MinMax).View=1;
Plugin(MinMax).OverTime=1;
Plugin(MinMax).Argument=0;
Plugin(MinMax).Run;
Save View [3] "tempMaxInterface.txt"; 


// EXTRACT TEMPERATURE FIELD AT PROBETIME
// Beware that ONELAB cannot always correctly guess the pulse endtime
// It is up to  the user to set PROBETIME 
Plugin(ExtractElements).MinVal=0;
Plugin(ExtractElements).MaxVal=0;
Plugin(ExtractElements).TimeStep = OL.eval(abs(OL.get(PostPro/2PROBETIME)/(OL.get(Parameters/Elmer/2TimeStep)*1e-3)));
Plugin(ExtractElements).Visible=1;
Plugin(ExtractElements).Dimension=2;
Plugin(ExtractElements).View=0;
Plugin(ExtractElements).Run;

//EVALUATE TEMPERATURE PROFILE AT INTERFACE
Plugin(CutPlane).A=0;
Plugin(CutPlane).B=-1;
Plugin(CutPlane).C=0;
Plugin(CutPlane).D=(OL.get(Parameters/Skin/3DERMIS)-1e-4)*1.e-3;
Plugin(CutPlane).ExtractVolume=0;
Plugin(CutPlane).RecurLevel=4;
Plugin(CutPlane).TargetError=0;
Plugin(CutPlane).View=4;
Plugin(CutPlane).Run;

//EVALUATE THE BOOLEAN FUNCTION: ABOVE THRESHOLD (1/0)
//MULTIPLY BY 2 PI X IN VIEW OF INTEGRATION
//MULTIPLY WITH 1E6 TO CONVERT IN MM2
Plugin(MathEval).Expression0= "Step(v0-OL.get(PostPro/5OVERTEMP))*2e6*pi*x";
Plugin(MathEval).TimeStep=-1;
Plugin(MathEval).View=5;
Plugin(MathEval).OtherTimeStep=-1;
Plugin(MathEval).OtherView=-1;
Plugin(MathEval).ForceInterpolation=0;
Plugin(MathEval).PhysicalRegion=-1;
Plugin(MathEval).Run;

//INTEGRATE
Plugin(Integrate).View=6;
Plugin(Integrate).Run; 

Save View [7] "Adelta.txt";
