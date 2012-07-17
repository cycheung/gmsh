OL.block
Khi.number(0.5, Parameters/,"Normalized hysteresis coefficient");
OL.endblock

General.RotationX = 40; 
General.RotationY = 0; 
General.RotationZ = 60; 

Plugin(MathEval).Expression0= "Atanh(Sqrt(v0*v0+v2*v2))*Sqrt(v0*v0+v2*v2)+0.5*Log(1-v0*v0-v2*v2)-(OL.get(Fields/Hx))*v0-(OL.get(Fields/Hz))*v2+OL.get(Khi)*Sqrt( (v0-(OL.get(Fields/Jpx)))^2+(v2-(OL.get(Fields/Jpz)))^2 )";
Plugin(MathEval).TimeStep=-1;
Plugin(MathEval).View=0;
Plugin(MathEval).OtherTimeStep=-1;
Plugin(MathEval).OtherView=-1;
Plugin(MathEval).ForceInterpolation=0;
Plugin(MathEval).PhysicalRegion=-1;
Plugin(MathEval).Run;
View[1].Name = "H=OL.get(Fields/Hx) (A/m), Jp=OL.get(Fields/Jpx) (T)";
View[1].NormalRaise = -1;
View[1].Axes = 0;
Save View [1] "Functional.pos" ;
Plugin(MinMax).View=1;
Plugin(MinMax).OverTime=0;
Plugin(MinMax).Argument=1;
Plugin(MinMax).Run;
Save View [2] "minimum.txt" ;
View[0].Visible=0;
View[2].Visible=0;
View[3].Visible=0;
View[1].ShowScale = 1; 
View[2].ShowScale = 0; 
View[3].ShowScale = 0; 



