l = 0.2;
r = 0.8;
R = l + r;

z = l / 2;

d  = 75;
cl = 2 * R / d;

Point(0) = {0, 0, -z, cl};

Point(1) = {-r, 0, -z, cl};
Point(2) = {-R, 0, -z, cl};
Point(3) = {+r, 0, -z, cl};
Point(4) = {+R, 0, -z, cl};

Circle(1) = {2, 0, 4};
Circle(2) = {1, 0, 3};

Line(3) = {2, 1};
Line(4) = {3, 4};

Line Loop(1) = {3, 2, 4, -1};
Plane Surface(1) = {1};

Extrude{0, 0, l}{
  Surface{1};
}

Physical Surface(5) = {13};
Physical Surface(6) = {21};
Physical Volume(7)  = {1};
