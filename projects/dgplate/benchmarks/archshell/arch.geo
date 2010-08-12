// Arch bending example

// Size
R = 10.;
L = 1.;
Lc1 = 0.01;

// Point
Point(1) = { -R  ,  0.  ,  0.  };
Point(2) = { -R  ,  0.  ,  L   };
Point(3) = {  0. ,  R   ,  L   };
Point(4) = {  0. ,  R   ,  0.  };
Point(5) = {  0. ,  0.  ,  0.  };
Point(6) = {  0. ,  0.  ,  L   };

// Arc and Line
Line(1) = {1,2};
Circle(2)  = {2,6,3};
Line(3) = {3,4};
Circle(4)  = {4,5,1}; 

// Surface definition
Line Loop(5) = {1,2,3,4};
Ruled Surface(6) = {5};

// Physical
Physical Surface(99) = {6};
Physical Line(11) = {1};
Physical Line(13) = {3};
Physical Point(333) = {3};
Physical Point(444) = {4};

// mesh
Transfinite Line {1,3} = 5 Using Progression 1;
Transfinite Line {2,4} = 41 Using Progression 1;
Transfinite Surface {6};
Recombine Surface {6};
Mesh.Smoothing = 100;
