#!/usr/bin/env python
#coding=utf-8

import numpy, pylab, math
import OnelabClient as onelab
OL = onelab.client()
print('\nStarting METAMODEL - Action = %s' %(OL.getString('python/Action')))

# parameters declaration
# the given value is a default value
dtime = OL.declareNumber('SRM/dtime', value=0.01, label='Time step [s]')
nstep = OL.declareNumber('SRM/nstep', value=150, label='Nb of time steps')
x = OL.declareNumber('IO/1POSITION', value=30, label='Rotor Position [deg]')
IR = OL.declareNumber('IO/4IR', value=1, label='Current phase R [A]') 
IS = OL.declareNumber('IO/5IS', value=0, label='Current phase R [A]')

# a depending variable
IT = OL.declareNumber('IO/6IT', value=1.-IR-IS, label='Current phase R [A]', readOnly=1)

# declarations with no value yields a readOnly parameter (highlighted in the Gui)
# This way, output quantities and readOnly parameters can be pre-declared
# and given a range, a label, etc
OL.declareNumber('IO/2MST', label='Torque (MST) [Nm]')
OL.declareNumber('IO/3VWP', label='Torque (VWP) [Nm]')

if OL.action == 'check' :
   OL.run('srm', 'srm.py ', '')  # check the FE model
   exit(0)

for i in range(15) :
   x = 3.0*i
   # a setNumber command overrules the value on server
   OL.setNumber('IO/1POSITION', value=x);
   OL.run('srm', 'srm.py ', '')  # run FE model

   mst = OL.getNumber('IO/2MST') # obtain value from Onelab server
   vwp = OL.getNumber('IO/3VWP')

   print "Angular position=%f Torque=%g %g" %(x,mst,vwp)







