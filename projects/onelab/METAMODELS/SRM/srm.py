import numpy
import onelab
import pylab
import math

rpm = 2.*math.pi/60.  # rpm to rad/s conversion
inertia = 8e-3;
dtime = 0.01; #seconds
nstep = 150;

#p = 3 number of stator pole pairs
#floor(x/30)%3 = 0 if x \in [ 0, 30 [ + k \pi/(2p)
#              = 1 if x \in [30, 60 [ + k \pi/(2p)
#              = 2 if x \in [60, 90 [ + k \pi/(2p)
def Current(angle):
   if (math.floor(angle/30.)%3) == 1:
      return 1
   else:
      return 0

def Cast(x,a,b):
   period = b-a
   if period > 0:
      x -= math.floor(x/period)*period
   return x

def Target(t):
   if(t<1):
      return 200
   else:
      return 100

def Load(t):
   if(t>0.5):
      return 0.15
   else:
      return 0.

times = [0.]  # s
angles = [0.] # deg
speeds = [0]  # rpm
msts = [0.]   # Nm
vwps = [0.]   # Nm

onelab.modelName("srm")

#pylab.clf()
pylab.ion()

pylab.subplot(211)
pylab.grid(True)
pylab.xlabel('Time [s]')
pylab.ylabel('Speed [rpm]')
lines, = pylab.plot(times, speeds)

pylab.subplot(212)
pylab.grid(True)
pylab.xlabel('Angle [deg]')
pylab.ylabel('Torque [Nm]')
lines2, = pylab.plot(angles, msts, '.-', markerfacecolor='red')

i=0
angle = angles[i]
speed = speeds[i]
time = i*dtime
while i<nstep:
   onelab.setNumber("IO/1POSITION", angle)
   alpha = 29. * math.tanh((Target(time)-speed)/10)
   IR = Current(angle-alpha)
   IS = Current(angle-alpha+30.)
   IT = Current(angle-alpha+60.)
   onelab.setNumber("IO/4IR", IR)
   onelab.setNumber("IO/5IS", IS)
   onelab.setNumber("IO/6IT", IT)
   onelab.metamodel("compute");
   mst = onelab.getNumber("IO/2MST")
   vwp = onelab.getNumber("IO/3VWP")
   #
   speed += ((mst-Load(time))/inertia*dtime)/rpm # rpm
   angle += ((speed*rpm)*dtime)/math.pi*180. # degrees
   i += 1
   time += dtime
   #
   times.append(time)
   angles.append(angle) # (Cast(angle,0,360))
   speeds.append(speed)
   msts.append(mst) 
   vwps.append(vwp)
   #
   print ("Currents=(%g, %g, %g) alpha=%g x=%g N=%g T=%g" %(IR,IS,IT,alpha,angle,speed,mst) )
   pylab.subplot(211)
   lines.set_xdata(times)
   lines.set_ydata(speeds)
   vmax = 1.05*max(math.fabs(v) for v in speeds)
   pylab.axis([0,nstep*dtime,-vmax,vmax])
   #
   pylab.subplot(212)
   lines2.set_xdata(angles)
   lines2.set_ydata(msts)
   tmax = 1.05*max(math.fabs(v) for v in msts)
   pylab.axis([0.,1080.,-tmax,tmax])
   #
   pylab.draw()






