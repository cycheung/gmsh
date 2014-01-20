r = 1;
R = 10;

d  = 10;
cl = r / d;

Point(0) = { 0,  0, 0, cl};

Point(1) = {+r,  0,  0, cl};
Point(2) = { 0, +r,  0, cl};
Point(3) = {-r,  0,  0, cl};
Point(4) = { 0, -r,  0, cl};
Point(5) = { 0,  0, +r, cl};
Point(6) = { 0,  0, -r, cl};

Circle(1) = {1, 0, 2};
Circle(2) = {2, 0, 3};
Circle(3) = {3, 0, 4};
Circle(4) = {4, 0, 1};

Circle(5) = {1, 0, 5};
Circle(6) = {5, 0, 3};
Circle(7) = {3, 0, 6};
Circle(8) = {6, 0, 1};

Circle(9)  = {2, 0, 5};
Circle(10) = {5, 0, 4};
Circle(11) = {4, 0, 6};
Circle(12) = {6, 0, 2};

Point(111) = {+R,  0,  0, cl};
Point(112) = { 0, +R,  0, cl};
Point(113) = {-R,  0,  0, cl};
Point(114) = { 0, -R,  0, cl};
Point(115) = { 0,  0, +R, cl};
Point(116) = { 0,  0, -R, cl};

Circle(111) = {111, 0, 112};
Circle(112) = {112, 0, 113};
Circle(113) = {113, 0, 114};
Circle(114) = {114, 0, 111};

Circle(115) = {111, 0, 115};
Circle(116) = {115, 0, 113};
Circle(117) = {113, 0, 116};
Circle(118) = {116, 0, 111};

Circle(119)  = {112, 0, 115};
Circle(1110) = {115, 0, 114};
Circle(1111) = {114, 0, 116};
Circle(1112) = {116, 0, 112};

Line Loop(1) = { 1,  9, -5};
Line Loop(2) = { 9,  6, -2};
Line Loop(3) = { 2,  7, 12};
Line Loop(4) = {12, -1, -8};
Line Loop(5) = {4,   5, 10};
Line Loop(6) = {10, -3, -6};
Line Loop(7) = {11,  8, -4};
Line Loop(8) = {7, -11, -3};

Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};
Ruled Surface(5) = {5};
Ruled Surface(6) = {6};
Ruled Surface(7) = {7};
Ruled Surface(8) = {8};

Line Loop(111) = { 111,  119, -115};
Line Loop(112) = { 119,  116, -112};
Line Loop(113) = { 112,  117, 1112};
Line Loop(114) = {1112, -111, -118};
Line Loop(115) = {114,   115, 1110};
Line Loop(116) = {1110, -113, -116};
Line Loop(117) = {1111,  118, -114};
Line Loop(118) = {117, -1111, -113};

Ruled Surface(111) = {111};
Ruled Surface(112) = {112};
Ruled Surface(113) = {113};
Ruled Surface(114) = {114};
Ruled Surface(115) = {115};
Ruled Surface(116) = {116};
Ruled Surface(117) = {117};
Ruled Surface(118) = {118};

Surface Loop(1) = {118, 113, 112, 111, 114, 117, 115, 116};
Surface Loop(2) = {6, 5, 7, 8, 3, 2, 1, 4};

Volume(1) = {1, 2};


Physical Surface(5) = {  1,   2,   3,   4,   5,   6,   7,   8};
Physical Surface(6) = {111, 112, 113, 114, 115, 116, 117, 118};
Physical Volume(7) = {1};
