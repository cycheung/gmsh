// Test case a SCB with a vertical load at its free extremity
// Size
x=0.05;
y=0.00625;

// Characteristic length
Lc1=0.001;

// definition of points
Point(1) = { 0.0 , 0.0 , 0.0 , Lc1};
Point(2) = {  x  , 0.0 , 0.0 , Lc1};
Point(3) = {  x  ,  y  , 0.0 , Lc1};
Point(4) = { 0.0 ,  y  , 0.0 , Lc1};
Point(5) = { x/2 , 0.0 , 0.0 , Lc1}; // Line at middle of beam to prescribed the displacement
Point(6) = { x/2 ,  y  , 0.0 , Lc1};

// Line between points
Line(1) = {1,5};
Line(2) = {5,2};
Line(3) = {2,3};
Line(4) = {3,6};
Line(5) = {6,4};
Line(6) = {4,1};
Line(7) = {5,6};

// Surface definition
Line Loop(8) = {1,7,5,6};
Line Loop(9) = {2,3,4,-7};
Plane Surface(10) = {8};
Plane Surface(11) = {9};
Physical Surface(99) = {11, 10};
Physical Line(61) = {6};
Physical Line(31) = {3};
Physical Line(71) = {7};
Transfinite Line {1, 2, 4, 5} = 5 Using Progression 1;
Transfinite Line {6, 7, 3} = 2 Using Progression 1;
Transfinite Surface {10};
Recombine Surface {10};
Transfinite Surface {11};
Recombine Surface {11};
Mesh.Smoothing = 100;
