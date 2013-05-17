cl  = 1;
msh = 75;

Point(1) = {-1, -1, 0, cl};
Point(2) = {+1, -1, 0, cl};
Point(3) = {+1, +1, 0, cl};
Point(4) = {-1, +1, 0, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1, 2, 3, 4} = msh Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1};

Physical Surface(7) = {1};
