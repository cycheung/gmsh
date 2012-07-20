cl = 0.2;

Point(1) = {0, 0, 0, cl};
Point(2) = {1, 0, 0, cl};
Point(3) = {0, 1, 0, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};

Line Loop(1) = {1, 2, 3};

Plane Surface(1) = {1};

Physical Line(5) = {1};
Physical Line(6) = {2};
Physical Surface(7) = {1};
