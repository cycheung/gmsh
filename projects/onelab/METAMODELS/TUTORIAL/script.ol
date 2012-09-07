OL.block
# place here ONELAB definitions
SATURATION.radioButton(1);
OL.endblock

# Read in a mesh
Merge "circle.msh";

# One could alternatively have the script generate the mesh with:
# Merge "circle.geo";
# Mesh 2;

#Create a(n empty) View[0] so that MathEval can be invoked
Plugin(NewView).Run; 

# Algebraic expression of the Functional
OL.iftrue(SATURATION)
Plugin(MathEval).Expression0= "Atanh(Sqrt(x*x+z*z))*Sqrt(x*x+z*z)+0.5*Log(1-x*x-z*z)-(OL.get(Fields/Hx))*x-(OL.get(Fields/Hz))*z+OL.get(Parameters/Khi)*Sqrt( (x-(OL.get(Fields/Jpx)))^2+(z-(OL.get(Fields/Jpz)))^2 )";
OL.else
Plugin(MathEval).Expression0= "(x*x+z*z)-(OL.get(Fields/Hx))*x-(OL.get(Fields/Hz))*z+OL.get(Parameters/Khi)*Sqrt( (x-(OL.get(Fields/Jpx)))^2+(z-(OL.get(Fields/Jpz)))^2 )";
OL.endif
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



