--[[
    Script to launch beam problem with a lua script
  ]]

-- creation of Solver
mysolver = DgC0PlateSolver(1000)
mysolver:readInputFile("beam.dat")
mysolver:createInterfaceElement()
mysolver:solve()
mysolver:buildDisplacementView("displacement")

