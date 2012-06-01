msh = 10;
cl = 1.5;

l = 1;

Point(1) = {+l, -l, 0, cl};
Point(2) = {+l, +l, 0, cl};
Point(3) = {-l, +l, 0, cl};
Point(4) = {-l, -l, 0, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Line {1, 2, 3, 4} = msh Using Progression 1;
Transfinite Surface {1};

Physical Line(5) = {3};
Physical Line(6) = {1};
Physical Surface(7) = {1};
