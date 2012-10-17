cl = 0.2;

Point(1) = {0, 0, 0, cl};
Point(2) = {1, 0, 0, cl};
Point(3) = {0, 1, 0, cl};
Point(4) = {0, 0, 1, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line(4) = {1, 4};
Line(5) = {2, 4};
Line(6) = {3, 4};

Line Loop(1) = {1, 2, 3};
Line Loop(2) = {2, 6, -5};
Line Loop(3) = {3, 4, -6};
Line Loop(4) = {1, 5, -4};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

Surface Loop(1) = {4, 1, 2, 3};
Volume(1)       = {1};

Physical Volume(7) = {1};
