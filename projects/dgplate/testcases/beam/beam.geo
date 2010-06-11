// Test case a SCB with a vertical load at its free extremity
// Size
x=0.1;
y=0.01;

// Characteristic length
Lc1=0.01;

// definition of points
Point(1) = { 0.0 , 0.0 , 0.0 , Lc1};
Point(2) = {  x  , 0.0 , 0.0 , Lc1};
Point(3) = {  x  ,  y  , 0.0 , Lc1};
Point(4) = { 0.0 ,  y  , 0.0 , Lc1};

// Line between points
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Surface definition
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

// Physical objects to applied BC and material
Physical Surface(99) = {6};
Physical Line(41) = {4};
Physical Line(31) = {3};
Physical Line(21) = {2};
Physical Line(11) = {1};
Physical Point(22) ={2};
Transfinite Line {2, 4} = 2 Using Progression 1;
Transfinite Line {1, 3} = 7 Using Progression 1;
Transfinite Surface {6};
Recombine Surface {6};
Mesh.Smoothing = 100;
