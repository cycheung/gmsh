import numpy
import onelab
import pylab
import math

nangle=180
mangle=180
dangle=mangle/nangle

angles = [0.]
msts = [0.]
vwps = [0.]

onelab.modelName("srm")
onelab.setNumber("IO/4IR", 1)
onelab.setNumber("IO/5IS", 0)
onelab.setNumber("IO/6IT", 0)

#l = plot(t1, f(t1), 'bo', t2, f(t2), 'k--', markerfacecolor='green')
pylab.ion()
pylab.grid(True)
pylab.xlabel('Angle [deg]')
pylab.ylabel('Torque [Nm]')
lines, = pylab.plot(angles, msts)

i=0
for i in range(0,nangle):
   angle = i*dangle
   onelab.setNumber("IO/1POSITION", angle)
   onelab.metamodel("compute");
   mst = onelab.getNumber("IO/2MST")
   vwp = onelab.getNumber("IO/3VWP")
   angles.append(angle)
   msts.append(mst) 
   vwps.append(vwp)
   print ("\n Position=%g Mst=%g Vwp=%g" %(angle,mst,vwp) )
   lines.set_xdata(angles)
   lines.set_ydata(msts)
   tmax = 1.05*max(math.fabs(v) for v in msts)
   pylab.axis([0,mangle,-tmax,tmax])
   pylab.draw()
   i=i+1

pylab.savefig("TA_new.pdf")


