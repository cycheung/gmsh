OL.block

LOADTYPE.number(0,2Load/2,"Type");
LOADTYPE.valueLabels(0, "None",
                     1, "Uniform",   
                     2, "Ramp", 
                     3, "Rectangle",
                     4, "Torsion",
                     5, "User defined");
FORMULA.string(,2Load/3,"Formula");

LOAD.number(1, 2Load/1,"Magnitude - P [N/cm2]");

D.number(0.1, 2Load/,"delta d [m]");
XX.number(0.5, 2Load/,"Position X [m]")
XX.range(0:OL.get(1Geometry/L):0.1);
XX.withinRange();

OL.if(OL.get(LOADTYPE) == 0)
FORMULA.setValue(0);
OL.endif
OL.if(OL.get(LOADTYPE) == 1)
FORMULA.setValue(P);
OL.endif
OL.if(OL.get(LOADTYPE) == 2)
FORMULA.setValue(P *tx(0)/L);
OL.endif
OL.if(OL.get(LOADTYPE) == 3)
FORMULA.setValue(if(tx(0)<(X-d)) {0} else {if(tx(0)<=(X+d)) {P/(2*d)} else {0}});
OL.endif
OL.if(OL.get(LOADTYPE) == 4)
FORMULA.setValue(P*tx(1)/B*10);
OL.endif

OL.if(OL.get(LOADTYPE) == 0)
LOAD.setVisible(0);
FORMULA.setVisible(0);
OL.else
LOAD.setVisible(1);
FORMULA.setVisible(1);
OL.endif

OL.if(OL.get(LOADTYPE) == 3)
D.setVisible(1);
XX.setVisible(1);
OL.else
D.setVisible(0);
XX.setVisible(0);
OL.endif

OL.if(OL.get(LOADTYPE) == 5)
FORMULA.setReadOnly(0);
OL.else
FORMULA.setReadOnly(1);
OL.endif

OL.endblock

Header
  Mesh DB "." "mesh"
End

Constants
End

Simulation
  Max Output Level = 30
  Coordinate System = Cartesian 3D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Steady State Min Iterations = 1
  Output Intervals = 1
  Post File = "elasticity.ep"
End

Solver 1
  Equation = "Displacement analysis"
  !Equation = "Elasticity solver"
  Variable = Displacement
  Variable DOFs = 3
  !Element = "p:2"
  Procedure = "StressSolve" "StressSolver"
  !Procedure = "ElasticSolve" "ElasticSolver"
  Linear System Solver = Iterative
  Linear System Iterative Method = CG
  Linear System Preconditioning = ILUT
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-6
  Linear System Residual Output = 20
  Nonlinear System Newton After Tolerance = 1.0e-3
  Nonlinear System Newton After Iterations = 20
  Nonlinear System Max Iterations = 1000
  Nonlinear System Convergence Tolerance = 1.0e-5
  Nonlinear System Relaxation Factor = 1.0
  Steady State Convergence Tolerance = 1.0e-4
  Calculate Loads = True
  Calculate stresses = True
End

Solver 2
   Exec solver = "After all"
   Equation = String "ResultOutput"
   Procedure = File "ResultOutputSolve" "ResultOutputSolver"
   Variable 1 = Stress
   Output File Name = String "beam.pos"
   Output Format = String "Gmsh"
   Gmsh Format = Logical true
End

Solver 3 !ElmerModelsManuel page 187
  Exec solver = "After all"
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Variable 1 = Displacement 2
  Operator 1 = min
  Operator 2 = max
  Save Coordinates(2,3) = 0 OL.eval(OL.get(1Geometry/A)/2) 0 OL.get(1Geometry/L) 0 0
  Filename = "beam.txt"
End

Equation 1
  Active Solvers(3) = 1 2 3
  Stress Analysis = True
End

Body OL.region(Volume)
  Equation = 1
  Material = 1
  Body Force = 1
End

Body Force 1
  Stress Bodyforce 1 = 0
  Stress Bodyforce 2 = Real MATC "-9.81 * OL.get(Material/DENSITY) * OL.get(Material/WEIGHT)"
  Stress Bodyforce 3 = 0
End 

Material 1
   Density = 1.0
   Youngs Modulus = Real OL.get(Material/YOUNG)e9
   Poisson Ratio = OL.get(Material/POISSON)
End

Boundary Condition 1
  Target Boundaries = OL.region(Clamping)
  Displacement 1 = 0
  Displacement 2 = 0
  Displacement 3 = 0
End

! variables from ONELAB
$P = OL.get(2Load/LOAD)*1e4
$L = OL.get(1Geometry/L)
$A = OL.get(1Geometry/A)
$B = OL.get(1Geometry/B)
$X = OL.get(2Load/XX)
$d = OL.get(2Load/D)

Boundary Condition 2
  Target Boundaries = OL.region(LoadSurf)
  Force 1 = 0
  Force 2 = Variable Coordinate 1, Coordinate 3
            Real MATC "OL.get(Load/FORMULA)"
  Force 3 = 0
End

Boundary Condition 3
  Target Boundaries = OL.region(FreeEnd)
  Force 1 = 0
  Force 2 = 0
End

