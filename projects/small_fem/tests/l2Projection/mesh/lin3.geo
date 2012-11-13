msh = 9;
cl  = 1;

Point(1) = {-0.5, 0, 0, cl};
Point(2) = {+0.5, 0, 0, cl};

Line(1) = {1, 2};

Transfinite Line {1} = msh Using Progression 1;

Physical Line(7) = {1};
