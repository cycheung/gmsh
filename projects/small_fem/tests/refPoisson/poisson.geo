l = 1;

msh = 3;
cl  = 5;

Point(1) = {+l, -l, 0, cl};
Point(2) = {+l, +l, 0, cl};
Point(3) = {-l, +l, 0, cl};
Point(4) = {-l, -l, 0, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1)     = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1, 2, 3, 4} = msh Using Progression 1;
Transfinite Surface {1};

Physical Line(5) = {1};
Physical Line(6) = {2};
Physical Line(7) = {3};
Physical Line(8) = {4};

Physical Surface(9) = {1};
