OL.block
PLANESTR.radioButton(1,1Geometry/,"Plane stress (1:0)");
PLANESTR.setVisible(0);
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
  Variable = Displacement
  Variable DOFs = 3
  !Element = p:2
  Procedure = "StressSolve" "StressSolver"
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU1
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
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
OL.iftrue(1Geometry/PLANESTR)
  Plane Stress = True
OL.else
  Plane Stress = False
OL.endif
End

Body OL.region(Volume)
  Equation = 1
  Material = 1
  Body Force = 1
End

Body Force 1
  Stress Bodyforce 1 = 0
  Stress Bodyforce 2 = Real MATC "-9.81 * OL.get(Material/DENSITY) * OL.get(Material/WEIGHT)"
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

Boundary Condition 2
  Target Boundaries = OL.region(LoadSurf)
  Force 1 = Variable Coordinate 1
            Real MATC "OL.get(Loads/LOADX) * tx(0)/OL.get(1Geometry/L)"
  Force 2 = Variable Coordinate 1
            Real MATC "OL.get(Loads/LOADY) * tx(0)/OL.get(1Geometry/L)"
  Force 3 = 0
End

Boundary Condition 3
  Target Boundaries = OL.region(LoadSurf)
  Force 1 = 0
  Force 2 = Variable Coordinate 1
            Real MATC " if(tx(0) > 0.5) {7000} else {-30000} "
  Force 3 = 0
End

Boundary Condition 4
  Target Boundaries = OL.region(FreeEnd)
  Force 1 = 0
  Force 2 = 0
End
