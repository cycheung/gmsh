L  = 1;
cl = 0.05;

Point(0) = {0, 0, 0, cl};

Point(1) = {+L,  0, 0, cl};
Point(2) = { 0, +L, 0, cl};
Point(3) = {-L,  0, 0, cl};
Point(4) = { 0, -L, 0, cl};

Circle(1) = {1, 0, 2};
Circle(2) = {2, 0, 3};
Circle(3) = {3, 0, 4};
Circle(4) = {4, 0, 1};

Line Loop(1)     = {3, 4, 1, 2};
Plane Surface(1) = {1};

Physical Line(5)    = {1, 2, 3, 4};
Physical Surface(7) = {1};
