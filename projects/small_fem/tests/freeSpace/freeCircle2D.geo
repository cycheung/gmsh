r = 1;
R = 10;

d  = 20;
cl = r / d;

Point(0) = { 0,  0, 0, cl};

Point(1) = {+r,  0, 0, cl};
Point(2) = { 0, +r, 0, cl};
Point(3) = {-r,  0, 0, cl};
Point(4) = { 0, -r, 0, cl};

Circle(1) = {1, 0, 2};
Circle(2) = {2, 0, 3};
Circle(3) = {3, 0, 4};
Circle(4) = {4, 0, 1};

Point(11) = {+R,  0, 0, cl};
Point(12) = { 0, +R, 0, cl};
Point(13) = {-R,  0, 0, cl};
Point(14) = { 0, -R, 0, cl};

Circle(11) = {11, 0, 12};
Circle(12) = {12, 0, 13};
Circle(13) = {13, 0, 14};
Circle(14) = {14, 0, 11};

Line Loop(1)     = { 1,  2,  3,  4};
Line Loop(2)     = {11, 12, 13, 14};
Plane Surface(1) = { 2,  1};

Physical Line(5)    = { 1,  2,  3,  4};
Physical Line(6)    = {11, 12, 13, 14};
Physical Surface(7) = {1};
