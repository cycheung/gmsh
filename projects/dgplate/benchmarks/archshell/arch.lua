--[[
    Script for arch bending benchmark
  ]]

--[[
    Data
  ]]

--material law
lawnum = 1
E = 3.e7
nu = 0.3

-- geometry
h = 1.
meshfile ="archO3.msh"
-- integration
nsimp = 1 

-- solver
sol = 1 
beta1 = 10.
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
mysolver:prescribedDisplacement("Edge",11,1,0.)
mysolver:prescribedDisplacement("Edge",13,0,0.)
--mysolver:prescribedDisplacement("Edge",11,2,0.)
mysolver:prescribedDisplacement("Face",99,2,0.)
mysolver:AddThetaConstraint(13)
mysolver:AddThetaConstraint(11)
mysolver:prescribedForce("Edge",13,0.,-5.e4,0.)

-- archiving
mysolver:ArchivingNodalDisplacement(3,1)
mysolver:ArchivingNodalDisplacement(4,1)
-- solve
mysolver:solve()
