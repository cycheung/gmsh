#!/usr/bin/env python
#coding=utf-8

import numpy, math
import onelab
OL = onelab.client()
print('\nStarting METAMODEL - Action = %s' %(OL.getString('python/Action')))

# Script calling the metamodel "srm.py" as a blackbox,
# here in a simple loop over the angular position

if OL.action == 'check' :
   OL.run('srm', './srm.py ', '')  # check the FE model
   exit(0)

for i in range(15) :
   x = 3.0*i
   # a setNumber command overrules the value on server
   OL.setNumber('IO/1POSITION', value=x);
   OL.run('srm', './srm.py ', '')  # run FE model

   mst = OL.getNumber('IO/2MST') # obtain value from Onelab server
   vwp = OL.getNumber('IO/3VWP')

   print "Angular position=%f Torque=%g %g" %(x,mst,vwp)







