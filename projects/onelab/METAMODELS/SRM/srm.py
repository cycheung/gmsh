#!/usr/bin/env python
#coding=utf-8

import OnelabClient as onelab
OL = onelab.client()


modelName = 'srm'

OL.geometry(modelName + '.geo')
if OL.action == 'compute' :
 OL.mesh(modelName + '.msh')

#OL.run('sub_gmsh', 'gmsh ' + modelName + '.geo -', '-2')


# call getDP, mesh file is given explicitly
# because the default mesh filename is guessed from the name of the main python script
# I.e., Info    : sub_solv - Got mesh name from Onelab: 'simple.msh'
# check command is the second argument
# run command is second + third argument
OL.run('sub_solv', 'getdp ' + modelName + '.pro', '-mesh ' + modelName + '.msh -solve MagSta')
OL.run('sub_post', 'getdp ' + modelName + '.pro', '-mesh ' + modelName + '.msh -pos Torque')



