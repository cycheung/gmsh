#-*-coding:Utf-8-*-
from math import pi
P = 1.e5
L = 1.
R = 10.
E = 3.e7
nu = 0.3
t = 1.
A = L*t
I = L*t*t*t/12.

dy = (P*L*R*R*R)/(2*E*I) * (pi/4- 2/pi) + (P*L*R*pi)/(8*E*A)
sigma = (P*R*L*t)/(2*pi*I)
print "la flèches est de : %.16f [m]"%dy
print "sigma max = %.16f [Pa]"%sigma
