
import numpy
import onelab

onelab.modelName("srm")

nangle = 10
dangle = 180/nangle
angle = numpy.zeros([nangle])
mst = numpy.zeros([nangle])
vwp = numpy.zeros([nangle])

onelab.setNumber("IO/4IR", 1)
onelab.setNumber("IO/5IS", 0)
onelab.setNumber("IO/6IT", 0)

filename=open('TorqueAngle.txt','w')
filename.write("angle mst vwp")

for i in range(0, nangle):
   angle = i*dangle
   print("The angle is %g (deg)" %(angle))
   onelab.setNumber("IO/1POSITION", angle[i])
   #onelab.metamodel("compute");
   mst[i] = onelab.getNumber("IO/2MST")
   vwp[i] = onelab.getNumber("IO/3VWP")
   print ("\n Position=%g Mst=%g Vwp=%g" %(angle[i],mst[i],vwp[i]) )
   i=i+1



