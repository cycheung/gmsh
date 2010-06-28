--[[
    Script to launch beam problem with a lua script
  ]]
--[[
    data 
  ]]
-- material law
lawnum = 1 -- unique number of the law
E = 100.e9 -- Young's modulus
nu = 0.   -- Poisson's ratio

-- geometry
h = 0.01  -- thickness
meshfile="beam50.msh" -- name of mesh file
-- integration
nsimp = 3 -- number of Simpson's points (odd)

-- solver
sol = 1 --Gmm=0 (default) Taucs=1 PETsc=2
beta1 = 10. -- value of stabilization parameter
beta2 = 10.
beta3 = 10.
soltype = 1 -- StaticLinear=0 (default) StaticNonLinear=1
nstep = 1   -- number of step (used only if soltype=1)
ftime =1.   -- Final time (used only if soltype=1)
tol=1.e-6   -- relative tolerance for NR scheme (used only if soltype=1)
nstepArch=1 -- Number of step between 2 archiving (used only if soltype=1)

--[[
    compute solution and BC (given directly to the solver
  ]]

-- creation of law
law1 = linearElasticLawPlaneStress(lawnum,E,nu)

-- creation of ElasticField
nfield =99 -- number of the field (physical number of surface)
fullDg = 0 --  formulation CgDg=0 fullDg =1
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
mysolver:prescribedDisplacement("Edge",41,0,0.)
mysolver:prescribedDisplacement("Edge",41,1,0.)
mysolver:prescribedDisplacement("Edge",41,2,0.)
--mysolver:prescribedDisplacement("Edge",21,0,0.)
--mysolver:prescribedDisplacement("Edge",21,1,0.)
--mysolver:prescribedDisplacement("Edge",21,2,0.)
--mysolver:prescribedDisplacement("Edge",21,1,0.4)
--mysolver:prescribedForce("Face",99,0.,0.,1000.)
mysolver:prescribedForce("Edge",21,0.,0.,-10000.)
mysolver:AddThetaConstraint(41)
mysolver:solve()
