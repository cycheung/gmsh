--[[
    Script for arch bending benchmark
  ]]

--[[
    Data
  ]]

--material law
lawnum = 1
E = 68.25e6
nu = 0.3

-- geometry
h = 0.04
meshfile = "sphereO2.msh"
-- integration
nsimp = 1 

-- solver
sol = 2 
beta1 = 1000.
beta2 = 10.
beta3 = 10.
soltype = 1 
nstep = 1
ftime = 1.
tol = 1.e-6
nstepArch = 1

--[[
    Compute solution and BC
  ]]

-- creation of material law
law1 = linearElasticLawPlaneStress(lawnum,E,nu)

-- creation of field
nfield = 99
fullDg = 1 
myfield1 = DGelasticField()
myfield1:tag(1000)
myfield1:thickness(h)
myfield1:lawnumber(lawnum)
myfield1:simpsonPoints(nsimp)
myfield1:formulation(fullDg)

-- creation of Solver
mysolver = DgC0PlateSolver(1000)
mysolver:readmsh(meshfile)
mysolver:AddElasticDomain(myfield1,nfield,2)
mysolver:AddLinearElasticLawPlaneStress(law1)
mysolver:setScheme(soltype)
mysolver:whichSolver(sol)
mysolver:SNLData(nstep,ftime,tol)
mysolver:stepBetweenArchiving(nstepArch)
mysolver:stabilityParameters(beta1,beta2,beta3)

-- BC
mysolver:prescribedDisplacement("Edge",12,0,0.)
mysolver:prescribedDisplacement("Edge",14,1,0.)
mysolver:prescribedDisplacement("Node",222,2,0.)

mysolver:prescribedForce("Node",222,0.,-20.,0.)
mysolver:prescribedForce("Node",111,20.,0.,0.)

mysolver:AddThetaConstraint(12)
mysolver:AddThetaConstraint(14)

-- archiving
mysolver:ArchivingNodalDisplacement(1,0)
mysolver:ArchivingNodalDisplacement(2,1)
-- solve
mysolver:solve()

