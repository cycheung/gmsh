msh = 10;
cl = 0.2;

Point(1) = {0, 0, 0, cl};
Point(2) = {1, 0, 0, cl};
Point(3) = {1, 1, 0, cl};
Point(4) = {0, 1, 0, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1)     = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1, 2, 3, 4} = msh Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1};

Extrude {0, 0, 1} {
  Surface{1};
  Layers{msh};
  Recombine;
}

Physical Line(5) = {3};
Physical Line(6) = {1};
Physical Volume(7) = {1};
