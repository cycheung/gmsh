
L = 0.004; 
H= 0.002;
R= 0.0015;

lc = 4*0.00002;

Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0., 0, lc};
Point(3) = {L,H, 0, lc};
Point(4) = {0, H, 0, lc};
Point(5) = {R, H, 0, lc};
Point(6) = {0, H+R, 0, lc};
 
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 5};
Line(4) = {5, 4};
Line(5) = {1, 4};
Circle(7) = {5, 4, 6};
Line(8) = {4, 6};
Line Loop(9) = {3, 4, -5, 1, 2};
Plane Surface(10) = {9};
Line Loop(11) = {7, -8, -4};
Plane Surface(12) = {11};

Physical Surface(1) = {10}; //peau
Physical Surface(2) = {12}; //verrue
Physical Line(11) = {1};
Physical Line(12) = {5};
Physical Line(13) = {2};
Physical Line(14) = {3};
Physical Line(15) = {7}; // cold
Physical Line(16) = {8};

