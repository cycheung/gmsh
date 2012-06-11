% in a name.xxx.ol file onelab parameter definition lines must be enclosed 
% between "OL.begin" and "OL.end" or start with "OL.line"
OL.begin
NumStep.number(10,Parameters/Elmer/1, "Time steps during laser appl."); 
TimeStep.number(OL.eval(OL.get(Parameters/Laser/APPLICTIME)/OL.get(NumStep)), Parameters/Elmer/2,"Time step [s]");
TimeStep.setReadOnly(1);
TimeEnd.number(OL.eval(3*OL.get(Parameters/Laser/APPLICTIME)),Parameters/Elmer/3,"Simulation end time [s]");
OL.end

%in the body of the file, onelab recognizes the following commands:
% OL.if, OL.if(n)true, OL.include, OL.eval, OL.get, OL.region

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
OL.if( OL.get(Parameters/Laser/LASERTYPE) == 3)
   Body Force = 2
OL.else
   Body Force = 3
OL.endif
End

Body 2 !Epidermis
  Equation = 1
  Material = OL.region(Epidermis)
  Initial Condition = 1
OL.if( OL.get(Parameters/Laser/LASERTYPE) == 3)
   Body Force = 1
OL.endif
End $

$teneurw  = OL.get(Parameters/Model/WCONTENT)
$pin      = OL.get(Parameters/Laser/LASERPOWER)
$r        = OL.get(Parameters/Model/BEAMRADIUS)/1000
$mua      = OL.get(Parameters/Laser/ABSORPTION)
$tlaser   = OL.get(Parameters/Laser/APPLICTIME)

$hp = (OL.get(Parameters/Model/DERMIS)+OL.get(Parameters/Model/SKINWIDTH))/1000
$ylaser =  hp

Body Force 1
Heat Source = Variable DensityBis, Coordinate 1, Coordinate 2, Time
Real MATC "if(tx(3)<=tlaser) {2*pin*(1-0.0078)/(pi*r*r)*mua*exp(-mua*(ylaser-tx(2))-2*tx(1)^2/(r*r))/tx(0)} else {0.0}"
End

Body Force 2
OL.if( OL.get(Parameters/Model/BIOHEAT) )
Heat Source = Variable DensityBis, Coordinate 1, Coordinate 2, Time, Temperature
Real MATC "if(tx(3)<=tlaser) {2*pin*(1-0.0078)/(pi*r*r)*mua*exp(-mua*(ylaser-tx(2))-2*(tx(1)/r)^2)/tx(0)+1452*(310-tx(4))} else {1452*(310-tx(4))}"
OL.else
Heat Source = Variable DensityBis, Coordinate 1, Coordinate 2, Time
Real MATC "if(tx(3)<=tlaser) {2*pin*(1-0.0078)/(pi*r*r)*mua*exp(-mua*(ylaser-tx(2))-2*tx(1)^2/(r*r))/tx(0)} else {0}"
OL.endif
End

Body Force 3
OL.if( OL.get(Parameters/Model/BIOHEAT) )
Heat Source = Variable DensityBis, Temperature
Real MATC "1452*(310-tx(1))/tx(0)"
OL.else
    Heat Source = Real 0.0
OL.endif
End

!*********** Equations ************
Equation 1
 Active Solvers(3) = 1 2 3
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
End
Solver 2
   Exec Solver = after saving
   Equation = String "ResultOutput"
   Procedure = File "ResultOutputSolve" "ResultOutputSolver"
   Output File Name = String "solution.pos"
   Output Format = String gmsh	
   Scalar Field 1 = String Temperature
   !Scalar Field 2 = String Teneur
   !Scalar Field 3 = String DensityBis
End
Solver 3 !ElmerModelsManuel page 187
   Exec solver = "After saving"
   Equation = "SaveScalars"
   Variable 1 = Temperature
   Variable 2 = Time
   Procedure = "SaveData" "SaveScalars"
   Save Coordinates(5,2) = 1e-6 OL.get(PostPro/ZSURF0) 1e-6 OL.get(PostPro/ZSURF1) 1e-6 OL.get(PostPro/ZSURF2) 1e-6 OL.get(PostPro/ZSURF3) 1e-6 OL.get(PostPro/ZSURF4)
   Filename = "temp.txt"
End


!*********** Materials ************
Material 1 !dermis
      OL.include( skinMaterial.ol )
End

Material 2 !epidermis
      OL.include( skinMaterial.ol )
End

Initial Condition 1
Temperature = Real OL.get(Parameters/Model/BODYTEMP)
Teneur = Variable Coordinate
OL.if( OL.get(Parameters/Model/SKINTYPE) == 1)
Real MATC "if((hp-tx(1))<0.00008){0.15/80*(hp-tx(1))*10^6+0.25} else {0.25+0.35/(1+exp(-0.25*((hp-tx(1))*10^6-80)))}"
OL.else
Real MATC "0.25+0.4/(1+exp(-0.25*((hp-tx(1))*10^6-15)))"
OL.endif
OL.if( OL.get(Parameters/Model/TENEUR) )
DensityBis = Variable Teneur
Real MATC "1000/(6.16/100*tx(0)+0.938)" ! kg/m3

HConductivity = Variable Teneur, DensityBis
Real MATC "tx(1)/1000*(0.454*tx(0)+0.174)" ! W/(mK) 

HCapacity = Variable Teneur
Real MATC "2500*tx(0)+1700"  ! J/(kg/K)
OL.else
DensityBis = Real MATC "1000/(6.16/100*teneurw+0.938)" !1022.45 !1048.88 ! 1200.0
OL.endif

End

!*********** Boundary conditions ************


Boundary Condition 1 ! "Zero flux on axis and bottom"
Target Boundaries(2) = OL.region(Axis) OL.region(Bottom)
Heat Flux BC = Logical true
Heat Flux Real = Real 0.0
End

Boundary Condition 2  ! "body temperature on far side"
Target Boundaries(1) = OL.region(Side)
Temperature = Real OL.get(Parameters/Model/BODYTEMP)
End

OL.if( OL.get(Parameters/Laser/LASERTYPE) == 1) 
   Boundary Condition 3  ! "laser spot"
   Target Boundaries(1) = OL.region(LaserSpot)
   Temperature = Variable Time 
   Real MATC "if(tx(0)<=tlaser) {OL.get(Parameters/Laser/LASERTEMP)} else {OL.get(Parameters/Model/BODYTEMP)}"
   End

   Boundary Condition 4 ! "convection or zero heat flux"
   OL.if( OL.get(Parameters/Model/CONVBC) )
   Target Boundaries(1) = OL.region(FreeSkin)
   Heat Transfer Coefficient = Real 75.0
   External Temperature = Real 292.15  ! 19 Celsius
   OL.else
   Heat Flux BC = Logical true
   Heat Flux Real = Real 0.0
   OL.endif
OL.endif

OL.if( OL.get(Parameters/Laser/LASERTYPE) == 2)
   Boundary Condition 3 ! "applied flux over whole skin surface"
   Target Boundaries(2) = OL.region(LaserSpot) OL.region(FreeSkin)
   Heat Flux BC = Logical true
   Heat Flux = Variable Coordinate 1, Coordinate 2, Time
   Real MATC "if(tx(2)<=tlaser) {2*pin*(1-0.0078)/(pi*r*r)*exp(-2*(tx(0)^2)/(r*r))} else {0.0}"
   End
OL.endif


OL.if( OL.get(Parameters/Laser/LASERTYPE) == 3)
   Boundary Condition 3 ! "applied zero flux over whole skin surface"
   Target Boundaries(2) = OL.region(LaserSpot) OL.region(FreeSkin)
   Heat Flux BC = Logical true
   Heat Flux = Real 0.0
OL.endif

