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

Point(11) = {-R  , 0.  , 0.  , Lc1};
Point(14) = {-R*sin, 0., R*cos,Lc1};

Point(22) = {0.,  -R  , 0.  , Lc1};
Point(23) = {0.,  -R*sin, R*cos,Lc1};

// Arc
Circle(1) = {1,5,2};
Circle(2) = {2,5,3};
Circle(3) = {3,6,4};
Circle(4) = {4,5,1};

Circle(11) = {11,5,2};
Circle(13) = {3,6,14};
Circle(14) = {14,5,11};

Circle(21) = {11,5,22};
Circle(23) = {23,6,14};
Circle(22) = {22,5,23};

Circle(31) = {1,5,22};
Circle(33) = {23,6,4};

// Surface definition
Line Loop(5) = {1,2,3,4};
Ruled Surface(6) = {5};

Line Loop(15) = {-11,-2,-13,-14};
Ruled Surface(16) = {15};

Line Loop(25) = {21,22,23,14};
Ruled Surface(26) = {25};

Line Loop(35) = {-31,-22,-33,-4};
Ruled Surface(36) = {35};


// Physical
Physical Surface(99) = {6,16,26,36};
Physical Line(11) = {1};
Physical Line(12) = {2};
Physical Line(13) = {3};
Physical Line(14) = {4};
Physical Point(111) = {1};
Physical Point(222) = {2};
Physical Point(1111) = {11};
Physical Point(2222) = {22};

// mesh
Transfinite Line {1,2,3,4,11,13,14,21,22,23,31,33} = 5;
Transfinite Surface {6,16,26,36};
Recombine Surface {6,16,26,36};
Mesh.Smoothing = 100;
