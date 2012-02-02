mm = 1.e-3;

DefineConstant[ H = {2*mm, Min 1*mm, Max 4*mm, Step 0.5*mm, Path "Parameters/1Geometry",ShortHelp "Skin thickness"} ];
DefineConstant[ L = {4*mm, Min 2*mm, Max 6*mm, Step 0.5*mm, Path "Parameters/1Geometry",ShortHelp "Model length"} ];
DefineConstant[ R = {2*mm, Min 1*mm, Max L-1*mm, Step 0.5*mm, Path "Parameters/1Geometry",ShortHelp "Radius"} ];

//DefineConstant[ lambda = {0.9, Min 0.7, Max 1, Step 0.05, Path "Parameters/1Geometry",ShortHelp "lambda"} ];
lambda=0.9;
DefineConstant[ Xloc = {lambda*R*Cos(Pi/4), Path "Parameters/1Geometry",ShortHelp "x coord of probepoint"} ];
DefineConstant[ Yloc = {lambda*R*Sin(Pi/4)+H, Path "Parameters/1Geometry",ShortHelp "y coord of probepoint"} ];

//DefineConstant[ Refine = {1, Min 0.01, Max 10, Path "Parameters/1Geometry",ShortHelp "Mesh Size"} ];
DefineConstant[ Refine = {1, Path "Parameters/1Geometry",ShortHelp "Mesh Size"} ];

lc=0.05*H*Refine;

Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0., 0, lc};
Point(3) = {L,H, 0, lc};
Point(4) = {0, H, 0, lc};
Point(5) = {R, H, 0, lc};
Point(6) = {0, H+R, 0, lc};
 
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 5};
Line(4) = {5, 4};
Line(5) = {1, 4};
Circle(7) = {5, 4, 6};
Line(8) = {4, 6};
Line Loop(9) = {3, 4, -5, 1, 2};
Plane Surface(10) = {9};
Line Loop(11) = {7, -8, -4};
Plane Surface(12) = {11};

Physical Surface(1) = {10}; //peau
Physical Surface(2) = {12}; //verrue
Physical Line(11) = {1};
Physical Line(12) = {5};
Physical Line(13) = {2};
Physical Line(14) = {3};
Physical Line(15) = {7}; // cold
Physical Line(16) = {8};

