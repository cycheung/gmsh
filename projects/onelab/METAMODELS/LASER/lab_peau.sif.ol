# in a name.xxx.ol file onelab parameter definition lines must be enclosed 
# between "OL.begin" and "OL.end" or start with "OL.line"
OL.block
NumStep.number(22,Parameters/Elmer/1, "Time steps during laser appl."); 
TimeStep.number(, Parameters/Elmer/2,"Time step [s]");
TimeStep.setValue(OL.eval(OL.get(Parameters/Laser/APPLICTIME)/OL.get(NumStep)));
TimeEnd.number(0.5,Parameters/Elmer/3,"Simulation end time [s]");
OL.endblock

#in the body of the file, onelab recognizes the following commands:
# OL.if, OL.if(n)true, OL.include, OL.eval, OL.get, OL.region

Header
  Mesh DB "." "mesh"
  echo on
End
Simulation
  Coordinate System =  Axi Symmetric 
  Simulation Type = Transient 
  Timestep Intervals (1) = OL.eval(floor(OL.get(TimeEnd)/OL.get(TimeStep)))
  Timestep sizes (1) = OL.get(TimeStep) 
  Timestepping Method = Implicit Euler
  ! BDF Order = 1
  Output Intervals = 1
End

!*********** Bodies ************
Body 1 !Dermis
  Equation = 1
  Material = OL.region(Dermis)
  Initial Condition = 1
  Body Force = 1
End

Body 2 !Epidermis
  Equation = 1
  Material = OL.region(Epidermis)
  Initial Condition = 1
  Body Force = 1
End

!*********** Equations ************
Equation 1
  Active Solvers(3) = 3 1 2
  !Convection = None
End
Solver 1
  Equation = "Heat Equation"
  Variable = String "Temperature"
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
  Exported Variable 1 = String "DensityBis"
  Exported Variable 1 DOFs = 1
  Exported Variable 2 = String "Teneur"
  Exported Variable 2 DOFs = 1
  Exported Variable 3 = String "Qvolume"
  Exported Variable 3 DOFs = 1
End
Solver 2
   Exec Solver = "After saving"
   Equation = String "ResultOutput"
   Procedure = File "ResultOutputSolve" "ResultOutputSolver"
   Output File Name = String "solution.pos"
   Output Format = String gmsh	
   Scalar Field 1 = String Temperature
   !Scalar Field 2 = String Teneur
   !Scalar Field 3 = String DensityBis
End
Solver 3 !ElmerModelsManuel page 187
  Exec solver = "Before timestep"
  Equation = "SaveScalars"
  Variable 1 = Temperature
  Variable 2 = Time
  Procedure = "SaveData" "SaveScalars"
  !Save Coordinates(5,2) = 1e-6 0.001549 1e-6 0.00150 1e-6 0.00145 
  ! 1e-6 0.0014 1e-6 0.00135 
  Save Coordinates(5,2) = 1e-6 OL.get(PostPro/ZSURF,choices.comp(0)) 1e-6 OL.get(PostPro/ZSURF,choices.comp(1)) 1e-6 OL.get(PostPro/ZSURF,choices.comp(2)) 1e-6 OL.get(PostPro/ZSURF,choices.comp(3)) 1e-6 OL.get(PostPro/ZSURF,choices.comp(4)) 
  Filename = "temp.txt"
  Target Variable 2 = String "Tsensor" 
End

!*********** Variables ************

$teneurw = OL.get(Parameters/Skin/WCONTENT)
$refl    = OL.get(Parameters/Skin/REFLECTIVITY)
$pin     = OL.get(Parameters/Laser/LASERPOWER)
$r       = OL.get(Parameters/Laser/BEAMRADIUS)/1000
$mua     = OL.get(Parameters/Laser/ABSORPTION)
$tlaser  = OL.get(Parameters/Laser/APPLICTIME)
$hp      = OL.get(PostPro/ZSURF)
$ylaser  = hp
$temp    = OL.get(Parameters/Laser/LASERTEMP)

!*********** Materials ************
Material 1 !dermis
# use a OL.include to avoid a duplicate section
OL.include( skinMaterial.ol )
End

Material 2 !epidermis
OL.include( skinMaterial.ol )
End

!*********** Initial condition ************

Initial Condition 1
  Temperature = Real OL.get(Parameters/Skin/BODYTEMP)
  Teneur = Variable Coordinate 1, Coordinate 2
OL.if( OL.get(Parameters/Skin/SKINTYPE) == 1) # hairy
  Real MATC "0.25+0.4/(1+exp(-0.25*((hp-tx(1))*10^6-15)))"
OL.else # hairless
  Real MATC "if((hp-tx(1))*10^6<80){0.15/80*(hp-tx(1))*10^6+0.25} else {0.25+0.35/(1+exp(-0.25*((hp-tx(1))*10^6-80)))}"
OL.endif

OL.if( OL.get(Parameters/Skin/TENEUR) )
  DensityBis = Variable Teneur
  Real MATC "1000/(0.0616*tx(0)+0.938)"
OL.else
  DensityBis = Real MATC "1000/(0.0616*teneurw+0.938)" 
OL.endif

  Qvolume = Variable DensityBis, Coordinate 1, Coordinate 2
OL.if( OL.get(Parameters/Laser/LASERSHAPE) == 1) # Gaussian
  Real MATC "(1-refl)*pin*2/(pi*r*r)*mua*exp(-mua*(ylaser-tx(2))-2*tx(1)^2/(r*r))/tx(0)"
OL.else # Flat-top
  Real MATC "if(tx(1)<r){(1-refl)*pin/(pi*r*r)*mua*exp(-mua*(ylaser-tx(2)))/tx(0)}else{0}"
OL.endif
End

!*********** Volume heat source ************
Body Force 1
OL.if( OL.get(Parameters/Laser/LASERTYPE) == 2)
  # imposed power density
  Heat Source = Variable Qvolume, Time
  Real MATC " if(tx(1)<=tlaser) {tx(0)} else {0} "
OL.endif
OL.if( OL.get(Parameters/Laser/LASERTYPE) == 3)
  # temperature controlled power density
  #Heat Source = Variable Qvolume, Time, Tsensor/Temperature
  Heat Source = Variable Qvolume, Time, Tsensor
  Real MATC " if(tx(1)<=tlaser) { if(tx(2)<temp){tx(0)} else{0} } else {0} "
OL.endif
End

!*********** Boundary conditions ************
Boundary Condition 1 ! "Zero flux on axis"
  Target Boundaries(1) = OL.region(Axis) 
  Heat Flux BC = Logical true
  Heat Flux Real = Real 0.0
End

Boundary Condition 2  ! "body temperature on external side and bottom"
  Target Boundaries(2) = OL.region(Side) OL.region(Bottom)
  Temperature = Real OL.get(Parameters/Skin/BODYTEMP)
End

OL.if( OL.get(Parameters/Laser/LASERTYPE) == 1) # imposed temperature
# In this case, 2 boundary conditions:
# imposed temperature on laser spot
Boundary Condition 3  ! 
  Target Boundaries(1) = OL.region(LaserSpot)
  Temperature = Variable Time 
  Real MATC "if(tx(0)<=tlaser) {OL.get(Parameters/Laser/LASERTEMP)} else {OL.get(Parameters/Skin/BODYTEMP)}"
End

# convection or zero heat on the free skin
Boundary Condition 4 ! 
  Target Boundaries(1) = OL.region(FreeSkin)
  OL.if( OL.get(Parameters/Skin/CONVBC) ) # convection
  Heat Transfer Coefficient = Real OL.get(Parameters/Skin/HCONV)
  External Temperature = Real OL.get(Parameters/Skin/TAMBIANT)
  OL.else # zero heat flux
  Heat Flux BC = Logical true
  Heat Flux Real = Real 0.0
  OL.endif
End

OL.else 
# In this case, power is injected through the volume source
# applied zero flux over whole skin surface
Boundary Condition 3 ! 
  Target Boundaries(2) = OL.region(LaserSpot) OL.region(FreeSkin)
  OL.if( OL.get(Parameters/Skin/CONVBC) ) # convection
  Heat Transfer Coefficient = Real OL.get(Parameters/Skin/HCONV)
  External Temperature = Real OL.get(Parameters/Skin/TAMBIANT)
  OL.else # zero heat flux
  Heat Flux BC = Logical true
  Heat Flux Real = Real 0.0
  OL.endif
End
OL.endif