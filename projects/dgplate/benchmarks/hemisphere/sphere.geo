// Pinched open hemisphere example

// Size
R = 10.;
Lc1 = 0.01;
// cos and sin value theta = 18°
sin = 0.3090169944;
cos = 0.9510565163;

// Point
Point(1) = { R  ,  0.  , 0.  , Lc1};
Point(2) = { 0. ,  R   , 0.  , Lc1};
Point(3) = { 0. , R*sin,R*cos, Lc1};
Point(4) = {R*sin, 0.  ,R*cos, Lc1};
Point(5) = { 0.  , 0.  , 0.  , Lc1};
Point(6) = { 0.  , 0.  ,R*cos, Lc1};

// Arc
Circle(1) = {1,5,2};
Circle(2) = {2,5,3};
Circle(3) = {3,6,4};
Circle(4) = {4,5,1};

// Surface definition
Line Loop(5) = {1,2,3,4};
Ruled Surface(6) = {5};

// Physical
Physical Surface(99) = {6};
Physical Line(11) = {1};
Physical Line(12) = {2};
Physical Line(13) = {3};
Physical Line(14) = {4};
Physical Point(111) = {1};
Physical Point(222) = {2};

// mesh
Transfinite Line {1,3} = 13;
Transfinite Line {2,4} = 13;
Transfinite Surface {6};
Recombine Surface {6};
Mesh.Smoothing = 100;
Line Loop(223) = {3, 4, 1, 2};
Ruled Surface(224) = {223};
