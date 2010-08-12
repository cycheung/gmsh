// Pinched Cylinder example

// Size
R = 0.3;
L = 0.6;
Lc1 = 0.01;

// Point
Point(1) = {  R  ,  0.  , 0.  , Lc1};
Point(2) = {  0. ,  R   , 0.  , Lc1};
Point(3) = {  0. ,  R   ,L/2. , Lc1};
Point(4) = {  R  ,  0.  ,L/2. , Lc1};
Point(5) = {  0. ,  0.  , 0.  , Lc1};
Point(6) = {  0. ,  0.  ,L/2. , Lc1};

// Arc and Line
Circle(1) = {1,5,2};
Line(2) = {2,3};
Circle(3) = {3,6,4};
Line(4) = {4,1};

// Surface definition
Line Loop(5) = {1,2,3,4};
Ruled Surface(6) = {5};

// Physical
Physical Surface(99) = {6};
Physical Line(11) = {1};
Physical Line(12) = {2};
Physical Line(13) = {3};
Physical Line(14) = {4};
Physical Point(333) = {3};
Physical Point(222) = {2};

// mesh
Transfinite Line {1,3} = 17;
Transfinite Line {2,4} = 17;
Transfinite Surface {6};
Recombine Surface {6};
Mesh.Smoothing = 100;  
  
