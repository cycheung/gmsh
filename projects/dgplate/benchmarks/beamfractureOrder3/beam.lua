--[[
    Script to launch beam problem with a lua script
  ]]
--[[
    data 
  ]]
-- material law
lawnum = 1
E = 71.e9 -- Young's modulus
nu = 0.   -- Poisson's ratio
Gc=8800.  -- fracture energy [J/m²]
sigmac=400.e6 -- fracture limit in tension [Pa]
beta =0.87 -- ratio KII/KI
mu =0.41 -- friction coefficient ??
-- geometry
h = 0.0025  -- thickness
meshfile="beamO3.msh" -- name of mesh file
-- integration
nsimp = 3 -- number of Simpson's points (odd)

-- solver
sol = 1 --Gmm=0 (default) Taucs=1 PETsc=2
beta1 =10. -- value of stabilization parameter
beta2 =10.
beta3 =10.
soltype = 1 -- StaticLinear=0 (default) StaticNonLinear=1
nstep =60   -- number of step (used only if soltype=1)
ftime =1.   -- Final time (used only if soltype=1)
tol=1.e-6   -- relative tolerance for NR scheme (used only if soltype=1)
nstepArch=1 -- Number of step between 2 archiving (used only if soltype=1)

--[[
    compute solution and BC (given directly to the solver
  ]]
-- creation of law
law1 = linearElasticLawPlaneStressWithFracture(lawnum,E,nu)
law1:setGc(Gc)
law1:setSigmac(sigmac)
law1:setBeta(beta)
law1:setMu(mu)

-- creation of ElasticField
nfield =99 -- number of the field (physical number of surface)
fullDg = 1 --  formulation CgDg=0 fullDg =1
myfield1 = dgLinearShellDomain()
myfield1:tag(1000)
myfield1:thickness(h)
myfield1:simpsonPoints(nsimp)
myfield1:formulation(fullDg)
myfield1:lawnumber(lawnum)
-- creation of Solver
mysolver = DgC0PlateSolver(1000)
mysolver:formulation(fullDg)
mysolver:readmsh(meshfile)
mysolver:addDgLinearElasticShellDomain(myfield1,nfield,2)
mysolver:AddLinearElasticLawPlaneStressWithFracture(law1)
mysolver:setScheme(soltype)
mysolver:whichSolver(sol)
mysolver:SNLData(nstep,ftime,tol)
mysolver:stepBetweenArchiving(nstepArch)
mysolver:stabilityParameters(beta1,beta2,beta3)
-- BC
mysolver:independentPrescribedDisplacement("Edge",31,0,0.00002)
mysolver:prescribedDisplacement("Edge",31,1,0.)
mysolver:prescribedDisplacement("Edge",31,2,0.0)
mysolver:prescribedDisplacement("Edge",61,0,0.)
mysolver:prescribedDisplacement("Edge",61,1,0.)
mysolver:prescribedDisplacement("Edge",61,2,0.)
mysolver:prescribedDisplacement("Edge",71,2,0.001)

mysolver:AddThetaConstraint(31)
mysolver:AddThetaConstraint(61)

mysolver:ArchivingForceOnPhysicalGroup(61,2)
mysolver:ArchivingForceOnPhysicalGroup(31,2)
mysolver:ArchivingForceOnPhysicalGroup(71,2)
mysolver:ArchivingNodalDisplacement(5,2)
mysolver:ArchivingNodalDisplacement(6,2)
mysolver:solve()

