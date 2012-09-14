cl = 2.75;
//cl = 1;
//cl = 0.75;

L = 2;
l = 0.75;

Point(0) = {0, 0, 0, cl};

Point(1) = {+L, -L, 0, cl};
Point(2) = {+L, +L, 0, cl};
Point(3) = {-L, +L, 0, cl};
Point(4) = {-L, -L, 0, cl};

Circle(1) = {1, 0, 2};
Circle(2) = {2, 0, 3};
Circle(3) = {3, 0, 4};
Circle(4) = {4, 0, 1};

Point(5) = {+l, -l, 0, cl};
Point(6) = {+l, +l, 0, cl};
Point(7) = {-l, +l, 0, cl};
Point(8) = {-l, -l, 0, cl};

Circle(5) = {5, 0, 6};
Circle(6) = {6, 0, 7};
Circle(7) = {7, 0, 8};
Circle(8) = {8, 0, 5};

Line Loop(9) = {3, 4, 1, 2};
Line Loop(10) = {7, 8, 5, 6};
Plane Surface(11) = {9, 10};

Physical Line(5) = {1, 2, 3, 4};
Physical Line(6) = {5, 6, 7, 8};
Physical Surface(7) = {11};
