#-*-coding:Utf-8-*-

# Analytical solution for plate

# Properties
L=10.
t=0.1
E=1e6
nu=0.3
P=200

# Compute D
D=E*t*t*t/(12*(1-nu*nu))

dzf=0.00561*P*L*L/D
dzt=0.007071*P*L*L/D
dzs=0.01160*P*L*L/D

# Solution
print "Analytical solution for All edges clamped : %f"%dzf
print "Analytical solution for two opposing edges clamed : %f"%dzt
print "Analytical solution for All edges simply supported : %f"%dzs

