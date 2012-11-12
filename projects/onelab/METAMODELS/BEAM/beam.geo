
DefineConstant[ L = {1, Path "1Geometry/", Label "Length"}];
DefineConstant[ A = {0.1, Path "1Geometry/", Label "Height"}];
DefineConstant[ B = {0.1, Path "1Geometry/", Label "Depth"}];
DefineConstant[ E = {0.05, Path "1Geometry/", Label "Characteristic element size", Visible 0}];
DefineConstant[CLAMPING = {1, Path "1Geometry/", 
			   Choices{1="One side", 2="Both sides"} }] ;

//DefineConstant [ RIV = {"RemoveInvisibleViews.geo", Path "Macros/", Label "Remove Invisible Views", Macro "Gmsh", Autocheck "0", Highlight "Turquoise"}];

NbLayX=20;
NbLayY=5;
NbLayZ=5;

If(E >= 0.001)
  NbLayX=Floor(L/E);
  NbLayY=2*Floor(B/E);
  NbLayZ=2*Floor(A/E);
EndIf

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
Coherence;

