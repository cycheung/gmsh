lc=0.03;
r =0.99;

Point(1) = { 0, 0, 0, lc};
Point(2) = { 0, r, 0, lc}; 
Line(1) = {1,2};

Point(3) = { r, 0, 0, lc};
Point(4) = { 0, 0, r, lc};
Point(5) = {-r, 0, 0, lc};
Point(6) = { 0, 0,-r, lc};

Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,6};
Circle(5) = {6,1,3};

Line Loop(6) = {2:5};
Plane Surface(7) = {6};
Physical Surface(1) = {7};
Physical Line(2) = {1};
