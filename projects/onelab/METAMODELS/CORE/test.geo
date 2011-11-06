/* --------------------------------------------------------------------------
    This is a sample gmsh geometry file
   --------------------------------------------------------------------------*/

/* defining some parameters */

e = 5.e-3 ;
d = 0.02 ;
h = 0.14 ;

DefineConstant[ l = {0.14, Min 0.05, Max 0.2, Step 0.01, Path "1Geometry",
                     ShortHelp "Core size"} ] ;
DefineConstant[ Val_Rint = {0.2, Min l, Max 50, Step 0.1, Path "1Geometry/1",
                            ShortHelp "Internal shell radius"}];
DefineConstant[ Val_Rext = {0.3, Min Val_Rint, Max 50, Step 0.1, Path "1Geometry/2",
                            ShortHelp "External shell radius"}];
ri = Val_Rint;
re = Val_Rext;

ha = 0.03 ;

p0 = d / 5 ;
p1 = e / 10 ;

pi = (re-ri)/8. ;

i = 0.001 ;
bob = 0.002 ;

xc = -l/2. ;

/* defining the points */

Point(1) = {xc+l/2,0,0,p0};
Point(2) = {xc+0,0,0,p0};
Point(3) = {xc+0,h/2,0,p0};
Point(4) = {xc+l,0,0,p1};
Point(5) = {xc+l,h/2,0,p0};
Point(6) = {xc+0,ha/2,0,p0};
Point(7) = {xc+d,ha/2,0,p0};
Point(8) = {xc+d,0,0,p0};
Point(9) = {xc+l-d,0,0,p1};
Point(10) = {xc+l-d,h/2-d,0,p0};
Point(11) = {xc+d,h/2-d,0,p0};
Point(12) = {xc+l,e/2,0,p1};
Point(13) = {xc+l-d,e/2,0,p1};
Point(14) = {xc+d+i,0,0,p0};
Point(15) = {xc+d+i,ha/2,0,p0};
Point(16) = {xc+d+i+bob,ha/2,0,p0};
Point(17) = {xc+d+i+bob,0,0,p0};
Point(18) = {xc+d+2*i+bob,0,0,p0};
Point(19) = {xc+d+2*i+bob,ha/2,0,p0};
Point(20) = {xc+d+2*i+2*bob,ha/2,0,p0};
Point(21) = {xc+d+2*i+2*bob,0,0,p0};
Point(22) = {xc-2*i-2*bob,0,0,p0};
Point(23) = {xc-2*i-2*bob,ha/2,0,p0};
Point(24) = {xc-2*i-bob,ha/2,0,p0};
Point(25) = {xc-2*i-bob,0,0,p0};
Point(26) = {xc-i-bob,0,0,p0};
Point(27) = {xc-i-bob,ha/2,0,p0};
Point(28) = {xc-i,ha/2,0,p0};
Point(29) = {xc-i,0,0,p0};
Point(30) = {xc+ri+l/2,0,0,pi};
Point(31) = {xc+re+l/2,0,0,pi};
Point(32) = {xc+l/2,ri,0,pi};
Point(33) = {xc+l/2,re,0,pi};
Point(34) = {xc+l/2-re,0,0,pi};
Point(35) = {xc+l/2-ri,0,0,pi};


/* defining the lines */

Line(1) = {22,23};
Line(2) = {23,24};
Line(3) = {24,25};
Line(4) = {25,22};
Line(5) = {26,29};
Line(6) = {29,28};
Line(7) = {28,27};
Line(8) = {27,26};
Line(9) = {2,6};
Line(10) = {8,7};
Line(11) = {7,6};
Line(12) = {2,8};
Line(13) = {14,17};
Line(14) = {17,16};
Line(15) = {16,15};
Line(16) = {15,14};
Line(17) = {18,21};
Line(18) = {21,20};
Line(19) = {20,19};
Line(20) = {19,18};
Line(21) = {7,11};
Line(22) = {6,3};
Line(23) = {3,5};
Line(24) = {10,11};
Line(25) = {10,13};
Line(26) = {13,9};
Line(27) = {4,12};
Line(28) = {12,5};
Line(29) = {21,1};
Line(30) = {1,9};
Line(31) = {9,4};
Line(32) = {12,13};
Line(33) = {2,29};
Line(34) = {2,8};
Line(35) = {25,26};
Line(36) = {34,35};
Line(37) = {35,22};
Line(38) = {31,30};
Line(39) = {33,32};
Line(40) = {4,30};
Circle(41) = {30,1,32};
Circle(42) = {32,1,35};
Circle(43) = {31,1,33};
Circle(44) = {33,1,34};
Line(51) = {8,14};
Line(52) = {17,18};

/* defining the surfaces */

Line Loop(45) = {39,43,-38,-41};
Plane Surface(46) = {45};

Line Loop(47) = {-42,-39,44,36};
Plane Surface(48) = {47};

Line Loop(49) = {-22,-23,28,-32,25,-24,21,-11};
Plane Surface(50) = {49};

Line Loop(53) = {37,1,2,3,35,-8,-7,-6,-33,9,22,23,-28,-27,40,41,42};
Plane Surface(54) = {53};

Line Loop(55) = {11,-9,34,10};
Plane Surface(56) = {55};

Line Loop(57) = {14,15,16,13};
Plane Surface(58) = {57};

Line Loop(59) = {18,19,20,17};
Plane Surface(60) = {59};

Line Loop(61) = {-1,-2,-3,-4};
Plane Surface(62) = {61};

Line Loop(63) = {6,7,8,5};
Plane Surface(64) = {63};

Line Loop(65) = {-16,-15,-14,52,-20,-19,-18,29,30,-26,-25,24,-21,-10,51};
Plane Surface(66) = {65};

Line Loop(67) = {31,27,32,26};
Plane Surface(68) = {67};


/* defining the physical entities (for which elements will be saved) */


Physical Surface(101) = {46,48} ;

Physical Surface(102) = {66,54} ;
Physical Surface(103) = {68} ;
Physical Surface(104) = {56} ;
Physical Surface(106) = {50} ;

Physical Surface(111) = {58} ;
Physical Surface(112) = {64} ;
Physical Surface(121) = {60} ;
Physical Surface(122) = {62} ;

Physical Line(1000) = {43,44} ;
Physical Line(1001) = {36,37,4,35,5,33,34,51,13,52,17,29,30,31,40,38};
