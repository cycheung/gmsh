OL.block
Tcold.number(77,Parameters/Elmer/2,"Applied temperature",50:100:10);
NumStep.number(50,Elmer/3,Number of time steps);
TimeStep.number(0.1,Elmer/3,Time step);
OL.endblock

Header
  Mesh DB "." "mesh"
End

Simulation
  Coordinate System =  Axi Symmetric
  Simulation Type = Transient 
  Timestep sizes = OL.get(TimeStep)
  Timestep Intervals = OL.get(NumStep)
  Timestepping Method = BDF
  BDF Order = 2
  Output Intervals = 1
  Solver Input File = "OL.get(Arguments/FileName).sif"
End


!*******************************
!*********** Bodies ************
!*******************************

Body 1
  Equation = 1
  Material = 1
  !Body Force = 1
  Initial Condition = 1
End

Body 2
  Equation = 1
  Material = 2
  Initial Condition = 1
End

!*******************************
!********* Equations ***********
!*******************************

Equation 1
  Name = "Heat Equation"
  Active Solvers(3) = 1 2 3
End

!*******************************
!*********** Solvers ***********
!*******************************


Solver 1
   Exec Solver = "Always"
   Equation = "Heat Equation"
   Variable = "Temperature"
   Variable Dofs = 1
   Linear System Solver = "Iterative"
   Linear System Iterative Method = "BiCGStab"
   Linear System Max Iterations = 350
   Linear System Convergence Tolerance = 1.0e-08
   Linear System Abort Not Converged = True
   Linear System Preconditioning = "ILU0"
   Linear System Residual Output = 1
   Steady State Convergence Tolerance = 1.0e-05
   Stabilize = True
   Nonlinear System Convergence Tolerance = 1.0e-05
   Nonlinear System Max Iterations = 1
   Nonlinear System Newton After Iterations = 3
   Nonlinear System Newton After Tolerance = 1.0e-02
   Nonlinear System Relaxation Factor = 1.0
   Adaptive Error Limit = 0.1
   Adaptive Max Change = 2
   Adaptive Coarsening = Logical True
End

Solver 2
   Exec solver = "After timestep"
   Equation = String "ResultOutput"
   Procedure = File "ResultOutputSolve" "ResultOutputSolver"
   Output File Name = String "solution.pos"
   Output Format = String "Gmsh"
   !Gmsh Format = Logical true	
   Scalar Field 1 = String Temperature	
End

Solver 3 !ElmerModelsManuel page 187
   Exec solver = "After saving"
   Equation = SaveScalars
   Variable 1 = Temperature
   Variable 2 = Time
   Procedure = "SaveData" "SaveScalars"
   Save Coordinates(1,2) = 0.0015 0.003
!These parameters were defined in the cryo.geo file
   Save Coordinates(1,2) = OL.get(Parameters/1Geometry/Xloc) OL.get(Parameters/1Geometry/Yloc)
   Filename = "tempevol.txt"
End

!*******************************
!***** Material Properties *****
!*******************************

Material 1
  Density = 1000
  Heat Conductivity = Variable Temperature
  Real
  18  0.627
  223 0.627
  248 0.62
  255 0.6
  260 0.55
  265 0.44
  270 0.27
  272 0.22
  273 0.209
  333 0.209
    End
  Heat Capacity =  Variable Temperature
  Real
  18    4180 
  248   4180
  261   5000
  265   10000
  268   20000
  269   80000
  270   44000
  270.5 20000
  271   4180
  333   4180
  End
End

Material 2
  Density = 1500
  Heat Conductivity = Variable Temperature
  Real
  18  0.627
  223 0.627
  248 0.62
  255 0.6
  260 0.55
  265 0.44
  270 0.27
  272 0.22
  273 0.209
  333 0.209
  End
  Heat Capacity =  Variable Temperature
  Real
  18    4180 
  248   4180
  261   5000
  265   10000
  268   20000
  269   80000
  270   44000
  270.5 20000
  271   4180
  333   4180
  End
End

!*******************************
!****** Intial Conditions ******
!*******************************

Initial Condition 1
 Temperature = Real 310.0
End

!*******************************
!**** Boundary Conditions ******
!*******************************

Boundary Condition 1
  Target Boundaries(1) = 15
  Heat Flux BC = Logical True
  Heat Transfer Coefficient = Real 5000. !Initial heat flux
  External Temperature = Real OL.get(Tcold)
End

Boundary Condition 2
  Target Boundaries(5) = 11 12 13 14 16 
  Heat Flux BC = Logical true
  Heat Flux Real = Real 0.0
End





