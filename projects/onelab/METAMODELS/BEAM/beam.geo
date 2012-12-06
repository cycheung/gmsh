
DefineConstant[ L = {1, Path "1Geometry/", Label "Length - L [m]"}];
DefineConstant[ A = {0.1, Path "1Geometry/", Label "Height - A [m]"}];
DefineConstant[ B = {0.1, Path "1Geometry/", Label "Width - B [m]"}];

DefineConstant[CLAMPING = {1, Path "1Geometry/", Label "Clamping",
			   Choices{1="One side", 2="Both sides"} }] ;
DefineConstant[STRUCTURED = {1, Path "1Geometry/", Label "Mesh",
			   Choices{0="Unstructured", 1="Structured"} }] ;

/*  == S T R U C T U R E D   M E S H ==  */

If(STRUCTURED == 1)
NbLayX=50;
NbLayY=7;
NbLayZ=7;
Printf("Number of layers = %g %g %g", NbLayX, NbLayY, NbLayZ);
lc=0.002; //dummy
Point(1) = {0,-A/2,-B/2, lc};
Point(2) = {0, A/2,-B/2, lc};
Point(3) = {0, A/2, B/2, lc};
Point(4) = {0,-A/2, B/2, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(10) = {5};

Transfinite Line {1, 3} = NbLayY;
Transfinite Line {2, 4} = NbLayZ;
Transfinite Surface {10} Right;
Recombine Surface {10};
Extrude Surface { 10, {L, 0, 0} } {Layers{NbLayX}; Recombine; } ;   

Physical Volume("Volume")={1};
Physical Surface("LoadSurf")={23};
Physical Surface("FreeEnd")={32};
If(CLAMPING == 1)
Physical Surface("Clamping")={10};
EndIf
If(CLAMPING == 2)
Physical Surface("Clamping")={10, 32};
EndIf
EndIf

/*  == U N S T R U C T U R E D   M E S H ==  */

If(STRUCTURED == 0)
lc1=(A+B)/15;
lc2=(A+B)/8;
lc3=lc2;
If(CLAMPING == 2)
lc3=lc1;
EndIf

Point( 1) = {0,-A/2,-B/2, lc1};
Point( 2) = {0, A/2,-B/2, lc1};
Point( 3) = {0, A/2, B/2, lc1};
Point( 4) = {0,-A/2, B/2, lc1};

Point( 5) = {L/2,-A/2,-B/2, lc2};
Point( 6) = {L/2, A/2,-B/2, lc2};
Point( 7) = {L/2, A/2, B/2, lc2};
Point( 8) = {L/2,-A/2, B/2, lc2};

Point( 9) = {L,-A/2,-B/2, lc3};
Point(10) = {L, A/2,-B/2, lc3};
Point(11) = {L, A/2, B/2, lc3};
Point(12) = {L,-A/2, B/2, lc3};

Line( 1) = { 1,  2};
Line( 2) = { 2,  3};
Line( 3) = { 3,  4};
Line( 4) = { 4,  1};

Line( 5) = { 1,  5};
Line( 6) = { 2,  6};
Line( 7) = { 3,  7};
Line( 8) = { 4,  8};

Line( 9) = { 5,  9};
Line(10) = { 6, 10};
Line(11) = { 7, 11};
Line(12) = { 8, 12};

Line(13) = { 9, 10};
Line(14) = {10, 11};
Line(15) = {11, 12};
Line(16) = {12,  9};

Line Loop(17) = {4, 1, 2, 3};
Plane Surface(18) = {17};
Line Loop(19) = {3, 8, 12, -15, -11, -7};
Plane Surface(20) = {19};
Line Loop(21) = {8, 12, 16, -9, -5, -4};
Plane Surface(22) = {21};
Line Loop(23) = {9, 13, -10, -6, -1, 5};
Plane Surface(24) = {23};
Line Loop(25) = {7, 11, -14, -10, -6, 2};
Plane Surface(26) = {25};
Line Loop(27) = {15, 16, 13, 14};
Plane Surface(28) = {27};

Surface Loop(29) = {18, 22, 20, 28, 24, 26};
Volume(30) = {29};


Physical Volume("Volume")={30};
Physical Surface("LoadSurf")={26};
Physical Surface("FreeEnd")={28};
If(CLAMPING == 1)
Physical Surface("Clamping")={18};
EndIf
If(CLAMPING == 2)
Physical Surface("Clamping")={18, 28};
EndIf
EndIf

EndIf
