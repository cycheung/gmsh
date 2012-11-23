
import numpy
import onelab
import string 
onelab.modelName("srm")

nangle = 45
dangle = 90/nangle

onelab.setNumber("IO/4IR", 1)
onelab.setNumber("IO/5IS", 0)
onelab.setNumber("IO/6IT", 0)

filename=open('TA.txt','w')
filename.write("angle[deg] Torque(mst) Torque(vwp) [Nm]\n")

for i in range(0, nangle):
   angle = i*dangle
   onelab.setNumber("IO/1POSITION", angle)
   onelab.metamodel("compute");
   mst = onelab.getNumber("IO/2MST")
   vwp = onelab.getNumber("IO/3VWP")
   filename.write(str(angle) + '\t' + str(mst) + '\t' + str(vwp) + '\n')
   i=i+1



