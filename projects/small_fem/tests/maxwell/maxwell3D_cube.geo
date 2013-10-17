msh = 3;
cl  = 1;
l   = 1;

L   = l / 2;

Point(1) = {+L, -L, -L, cl};
Point(2) = {+L, +L, -L, cl};
Point(3) = {-L, +L, -L, cl};
Point(4) = {-L, -L, -L, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1)     = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1, 2, 3, 4} = msh Using Progression 1;
Transfinite Surface {1};
Recombine Surface{1};

Extrude {0, 0, 2 * L} {
  Surface{1};
  Layers{msh - 1};
  Recombine;
}

Physical Surface(5) = {21};
Physical Surface(6) = {13};
Physical Volume(7)  = {1};
