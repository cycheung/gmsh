#!/usr/bin/env python
#coding=utf-8

import OnelabClient as onelab
OL = onelab.client()

# Script calculating the torque of a SRM
# as a function of phase currents and rotor angular position.

# declaration of the metamodel parameters
# the given value is a default value
dtime = OL.defineNumber('SRM/dtime', value=0.01, label='Time step [s]')
nstep = OL.defineNumber('SRM/nstep', value=150, label='Nb of time steps')
x = OL.defineNumber('IO/1POSITION', value=30, label='Rotor Position [deg]')
IR = OL.defineNumber('IO/4IR', value=1, label='Current phase R [A]') 
IS = OL.defineNumber('IO/5IS', value=0, label='Current phase R [A]')

# declaration of a depending variable (readOnly=1)
IT = OL.defineNumber('IO/6IT', value=-IR-IS, label='Current phase R [A]', readOnly=1)

# Output quantities can be declared with no value
# They are automatically made readOnly and highlighte in the Gui
# This declaration causes the value be reset by a "Check"
OL.defineNumber('IO/2MST', label='Torque (MST) [Nm]')
OL.defineNumber('IO/3VWP', label='Torque (VWP) [Nm]')

modelName = 'srm'

# 2 ways to call the client gmsh (my favourite is 2)

# 1. call gmsh as a client 
## OL.run('sub_gmsh', 'gmsh ' + modelName + '.geo -', '-2')

# 2. use OnelabClient.py specific functions
OL.openGeometry(modelName + '.geo')  # merge geo file
if OL.action == 'compute' :
  OL.mesh(modelName + '.msh') # send string "Mesh 3; Save filename;"

# call getDP, mesh file is given explicitly
# because the default mesh filename is guessed from the name of the main python script
# I.e., Info    : sub_solv - Got mesh name from Onelab: 'simple.msh'
# check command is the second argument
# run command is second + third argument
OL.run('sub_solv', 'getdp ' + modelName + '.pro', '-mesh ' + modelName + '.msh -solve MagSta')
OL.run('sub_post', 'getdp ' + modelName + '.pro', '-mesh ' + modelName + '.msh -pos Torque')



