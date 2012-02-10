
L = 0.004; 
H= 0.002;
lc = 4*0.00002;

Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0., 0, lc};
Point(3) = {L,H, 0, lc};
Point(4) = {0, H, 0, lc};
Point(5) = {H, H, 0, lc};
Point(6) = {0, L, 0, lc};
 
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 5};
Line(4) = {5, 4};
Line(5) = {1, 4};
Circle(6) = {5, 4, 6};
Line(7) = {4, 6};

Line Loop(1) = {1, 2, 3, 4,-5};
Plane Surface(1) = {1};
Line Loop(2) = {4, 7, -6};
Plane Surface(2) = {2};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {6};
Physical Line(5) = {7};
Physical Line(6) = {5};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
