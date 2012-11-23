// Switched Reluctance Motor Gmsh geometry file (2D)
// All dimensions in meters and rads

Include "srm_data.pro" ;

Mesh.Algorithm = 1 ;
Mesh.CharacteristicLengthFactor = 1. ;

// Some general settings
Lc = 0.00075 * 0.5 ;          // Base char length

// Stator data
th0ss=0.;       	// Angular pos. stator

Betas=32./180.*Pi;	// Opening Angle Stator Pole
dths=2.*Pi/Ns;		// Ang. shift between 2 poles

Rsout=0.06;		// Ext. diam.
Rsin=0.030;		// Int. diam.
YWs=0.012;		// Yoke width

// Airgap
AG=0.29e-3;

// Rotor data
angleR=0;
th0rs = angleR*Pi/180 ; // Angular pos. rotor

Betar = 32/180*Pi;	// Opening Angle Rotor Pole
dthr  = 2.*Pi/Nr;	// Ang. shift between 2 poles

Rrout=Rsin-AG;		// Ext. diam.
Rrin=0.02;		// Int. diam.
Rshaft=0.006;		// Shaft diam.

//------------------------------------------------------------
//------------------------------------------------------------

// Create origin
pAxe = newp ; Point(pAxe) = { 0. , 0. , 0., 3*Lc} ;

fac = 2 ; // Half model
// Create all stator points
Rsmid = Rsout-YWs ;
For N In {0:Ns/fac-1}
  th0s=N*dths+th0ss;
  p1sox=Rsout*Cos(-dths/2.+th0s); p1soy=Rsout*Sin(-dths/2.+th0s);
  p1six=Rsmid*Cos(-dths/2.+th0s); p1siy=Rsmid*Sin(-dths/2.+th0s);

  th4s=Asin(Rsin/Rsmid*Sin(Betas/2.));
  p2six=Rsmid*Cos(-th4s+th0s);    p2siy=Rsmid*Sin(-th4s+th0s);
  p3six=Rsin*Cos(-Betas/2.+th0s); p3siy=Rsin*Sin(-Betas/2.+th0s);
  p4six=Rsin*Cos( Betas/2.+th0s); p4siy=Rsin*Sin( Betas/2.+th0s);
  p5six=Rsmid*Cos(th4s+th0s);     p5siy=Rsmid*Sin(th4s+th0s);

  p1so[N]=newp; Point(p1so[N])={p1sox,p1soy,0.,8*Lc} ;
  p1si[N]=newp; Point(p1si[N])={p1six,p1siy,0.,8*Lc} ;
  p2si[N]=newp; Point(p2si[N])={p2six,p2siy,0.,8*Lc} ;
  p3si[N]=newp; Point(p3si[N])={p3six,p3siy,0.,1*Lc} ;
  p4si[N]=newp; Point(p4si[N])={p4six,p4siy,0.,1*Lc} ;
  p5si[N]=newp; Point(p5si[N])={p5six,p5siy,0.,8*Lc} ;
EndFor

N = Ns/fac ;
th0s=N*dths+th0ss;
p1sox=Rsout*Cos(-dths/2.+th0s);  p1soy=Rsout*Sin(-dths/2.+th0s);
p1six=Rsmid*Cos(-dths/2.+th0s);  p1siy=Rsmid*Sin(-dths/2.+th0s);
p1so[]+=newp; Point(newp)={p1sox,p1soy,0.,8*Lc} ;
p1si[]+=newp; Point(newp)={p1six,p1siy,0.,8*Lc} ;

// Create Stator Lines, arcs and regions
// outer stator surface
For N In {0:Ns/fac-1}
 arcso[]+=newl; Circle(newl)={p1so[N],pAxe,p1so[(N+1)%Ns]};
EndFor

// outer surface of N-th stator tooth
For N In {0:Ns/fac-1}
 clsi1[N]=newl; Circle(clsi1[N])={p1si[N],pAxe,p2si[N]};
 clsi2[N]=newl; Line(clsi2[N])={p2si[N],p3si[N]};
 clsi3[N]=newl; Circle(clsi3[N])={p3si[N],pAxe,p4si[N]};
 clsi4[N]=newl; Line(clsi4[N])={p4si[N],p5si[N]};
 clsi5[N]=newl; Circle(clsi5[N])={p5si[N],pAxe,p1si[(N+1)%Ns]};
EndFor

ss1=newl; Line(newl)={p1si[0], p1so[0]};
ss2=newl; Line(newl)={p1si[Ns/fac], p1so[Ns/fac]};
llStator=newll;
statorin[] = {clsi1[0]:clsi5[Ns/fac-1]} ;

llstator = newll ; Line Loop (llstator)={-ss1,clsi1[0]:clsi5[Ns/fac-1],ss2,-arcso[{2:0:-1}]};
Stator[] += news ; Plane Surface(Stator[0]) = {llstator} ;

// Create Coil regions
For N In {0:Ns/fac}
 th0s=N*dths+th0ss;
 p1cx=Rsin*Cos(-dths/2.+th0s);
 p1cy=Rsin*Sin(-dths/2.+th0s);
 p1c[N]=newp; Point(p1c[N])={p1cx,p1cy,0.,2*Lc} ;
EndFor
For N In {0:Ns/fac}
 clci1[N]=newl; Line(clci1[N])={p1si[N],p1c[N]};
EndFor
For N In {0:Ns/fac-1}
   clci2[N]=newl; Circle(clci2[N])={p1c[N],pAxe,p3si[N]};
   clci3[N]=newl; Circle(clci3[N])={p4si[N],pAxe,p1c[(N+1)%Ns]};
EndFor

For N In {0:Ns/fac-1}
 Coillln[N]=newll; Line Loop (Coillln[N])={clci1[N],clci2[N],-clsi2[N],-clsi1[N]};
 Coiln[N]=news ;   Plane Surface(Coiln[N])= {Coillln[N]};
 Coilllp[N]=newll; Line Loop (Coilllp[N])={-clsi4[N],clci3[N],-clci1[(N+1)%Ns],-clsi5[N]};
 Coilp[N]=news ;   Plane Surface(Coilp[N])= {Coilllp[N]};
EndFor


// Create all rotor points
For N In {0:Nr/fac-1}
  th0r = N*dthr+th0rs;
  p1rox=Rrin*Cos(-dthr/2.+th0r);   p1roy=Rrin*Sin(-dthr/2.+th0r);
  p6rox=Rrout*Cos(-dthr/2.+th0r);  p6roy=Rrout*Sin(-dthr/2.+th0r);

  th2r=Asin(Rrout/Rrin*Sin(Betar/2.));
  p2rox=Rrin*Cos(-th2r+th0r);      p2roy=Rrin*Sin(-th2r+th0r);
  p3rox=Rrout*Cos(-Betar/2.+th0r); p3roy=Rrout*Sin(-Betar/2.+th0r);
  p4rox=Rrout*Cos(Betar/2.+th0r);  p4roy=Rrout*Sin(Betar/2.+th0r);
  p5rox=Rrin*Cos(th2r+th0r);       p5roy=Rrin*Sin(th2r+th0r);
  p1rix=Rshaft*Cos(-dthr/2.+th0r); p1riy=Rshaft*Sin(-dthr/2.+th0r);

  p1ro[N]=newp; Point(p1ro[N])={p1rox,p1roy,0.,4*Lc};
  p2ro[N]=newp; Point(p2ro[N])={p2rox,p2roy,0.,6*Lc};
  p3ro[N]=newp; Point(p3ro[N])={p3rox,p3roy,0.,2*Lc/2};
  p4ro[N]=newp; Point(p4ro[N])={p4rox,p4roy,0.,2*Lc/2};
  p5ro[N]=newp; Point(p5ro[N])={p5rox,p5roy,0.,6*Lc};
  p6ro[N]=newp; Point(p6ro[N])={p6rox,p6roy,0.,6*Lc};
  p1ri[N]=newp; Point(p1ri[N])={p1rix,p1riy,0.,6*Lc};
EndFor

N = Nr/fac ;
th0r = N*dthr+th0rs;
p1rox=Rrin*Cos(-dthr/2.+th0r);   p1roy=Rrin*Sin(-dthr/2.+th0r);
p1rix=Rshaft*Cos(-dthr/2.+th0r); p1riy=Rshaft*Sin(-dthr/2.+th0r);
p6rox=Rrout*Cos(-dthr/2.+th0r);  p6roy=Rrout*Sin(-dthr/2.+th0r);
p1ro[N]=newp; Point(p1ro[N])={p1rox,p1roy,0.,4*Lc};
p1ri[N]=newp; Point(p1ri[N])={p1rix,p1riy,0.,6*Lc};
p6ro[N]=newp; Point(p6ro[N])={p6rox,p6roy,0.,6*Lc};

// Create Rotor Lines, arcs and regions
For N In {0:Nr/fac-1}
  arcri[N]=newl; Circle(arcri[N])={p1ri[N],pAxe,p1ri[(N+1)%Nr]};
EndFor
cutshaft[0] = newl; Line(newl)={p1ri[2],pAxe};
cutshaft[1] = newl; Line(newl)={pAxe,p1ri[0]};

For N In {0:Nr/fac-1}
  clro1[N]=newl; Circle(clro1[N])={p1ro[N],pAxe,p2ro[N]};
  clro2[N]=newl; Line(clro2[N])={p2ro[N],p3ro[N]};
  clro3[N]=newl; Circle(clro3[N])={p3ro[N],pAxe,p4ro[N]};
  clro4[N]=newl; Line(clro4[N])={p4ro[N],p5ro[N]};
  clro5[N]=newl; Circle(clro5[N])={p5ro[N],pAxe,p1ro[(N+1)%Nr]};
EndFor

rr1=newl; Line(newl)={p1ri[0],p1ro[0]};
rr2=newl; Line(newl)={p1ri[Nr/fac],p1ro[Nr/fac]};
llshaft = newll ; Line Loop (llshaft) = {arcri[],cutshaft[]};
Shaft[] += news ; Plane Surface(Shaft[0]) = {llshaft} ;

rotorout[]= clro1[0]:clro5[Nr/fac-1] ;
llrotor = newll ; Line Loop (llrotor) = {rotorout[],-rr2, -arcri[{1:0}],rr1};
Rotor[] += news ; Plane Surface(Rotor[0]) = {llrotor} ;
rr1_=newl; Line(newl)={p1ro[0],p6ro[0]};
rr2_=newl; Line(newl)={p1ro[2],p6ro[2]};

For N In {0:Nr/fac-1}
  clrr2[N]=newl; Circle(clrr2[N])={p6ro[N],pAxe,p3ro[N]};
  clrr3[N]=newl; Circle(clrr3[N])={p4ro[N],pAxe,p6ro[(N+1)%Ns]};
EndFor

Line Loop(newll) = {rotorout[{0,1}],-clrr2[0], -rr1_};
AirgapRotorIn[]+=news; Plane Surface(news) = {newll-1};
Line Loop(newll) = {rotorout[{3:6}], -clrr2[1],-clrr3[0]};
AirgapRotorIn[]+=news; Plane Surface(news) = {newll-1};
Line Loop(newll) = {rotorout[{8,9}], rr2_, -clrr3[1]};
AirgapRotorIn[]+=news; Plane Surface(news) = {newll-1};


//============================================================
// moving band - from Stator
Rag  = Rsin -AG/3 ; // radious from inner stator radious
Ragr = Rrout+AG/3 ; // radious from outer rotor radious

// Lines limiting the stator and coils (outer airgap)
For N In {0:Ns/fac-1}
  airgco[3*N]=clci2[N] ;
  airgco[3*N+1]=clsi3[N] ;
  airgco[3*N+2]=clci3[N] ;
EndFor

For N In {0:Ns/fac-1}
  th0s=N*dths+th0ss;
  p1mbsx=Rag*Cos(-dths/2.+th0s);  p1mbsy=Rag*Sin(-dths/2.+th0s);
  p2mbsx=Rag*Cos(-Betas/2.+th0s); p2mbsy=Rag*Sin(-Betas/2.+th0s);
  p3mbsx=Rag*Cos( Betas/2.+th0s); p3mbsy=Rag*Sin( Betas/2.+th0s);
  pmbs[]+=newp; Point(newp)={p1mbsx,p1mbsy,0.,1*Lc} ;
  pmbs[]+=newp; Point(newp)={p2mbsx,p2mbsy,0.,1*Lc} ;
  pmbs[]+=newp; Point(newp)={p3mbsx,p3mbsy,0.,1*Lc} ;
  If(N==Ns/fac-1)
    p4mbsx=Rag*Cos( dths/2.+th0s);  p4mbsy=Rag*Sin( dths/2.+th0s);
    pmbs[]+=newp; Point(newp)={p4mbsx,p4mbsy,0.,1*Lc} ;
  EndIf
EndFor
For N In {0:#pmbs[]-2}
  clmbs[]+=newl ; Circle(newl)={pmbs[N], pAxe, pmbs[N+1]} ;
EndFor
For N In {0:#clmbs[]-1}
  clmbs[]+= Rotate {{0, 0, 1}, {0, 0, 0}, Pi} { Duplicata { Line{clmbs[N]}; } };
EndFor

lmbs[]+=newl; Line(newl)={pmbs[0],p1c[0]};
lmbs[]+=newl; Line(newl)={p1c[Ns/fac],pmbs[#pmbs[]-1]};

Nl = #clmbs[] ;
llmbs=newll; Line Loop(llmbs)={airgco[],lmbs[1],-clmbs[{Nl/2-1:0:-1}],lmbs[0]};
surfmbstator[]+=news ; Plane Surface(surfmbstator[0]) = {llmbs};


// moving band - from Rotor
For N In {0:Nr-1}
  th0r = N*dthr+th0rs;
  p1mbrx= Ragr*Cos(-dthr/2.+th0r);  p1mbry=Ragr*Sin(-dthr/2.+th0r);
  p2mbrx= Ragr*Cos(-Betar/2.+th0r); p2mbry=Ragr*Sin(-Betar/2.+th0r);
  p3mbrx= Ragr*Cos( Betar/2.+th0r); p3mbry=Ragr*Sin( Betar/2.+th0r);
  pmbr[]+=newp; Point(newp)={p1mbrx,p1mbry,0.,2*Lc};//4
  pmbr[]+=newp; Point(newp)={p2mbrx,p2mbry,0.,2*Lc/2};
  pmbr[]+=newp; Point(newp)={p3mbrx,p3mbry,0.,2*Lc/2};
EndFor

For N In {0:#pmbr[]-1}
  clmbr[]+=newl ; Circle(newl)={pmbr[N], pAxe, pmbr[{(N<#pmbr[]-1)?N+1:0}]} ;
EndFor
lmbr[]+=newl; Line(newl)={pmbr[0],p6ro[0]};
lmbr[]+=newl; Line(newl)={p6ro[2],pmbr[6]};

bndmbrotor[] = {clrr3[1],rotorout[7],clrr2[1],clrr3[0],rotorout[2],clrr2[0]};
llmbr=newll; Line Loop(llmbr)={-lmbr[0],clmbr[{0:#clmbr[]/2-1}],-lmbr[1], -bndmbrotor[]};
surfmbrotor[]+=news ; Plane Surface(surfmbrotor[0]) = {llmbr};


llmbin = newll ; Line Loop(llmbin)= {clmbr[]};
llmbout = newll ; Line Loop(llmbout)= {clmbs[]};
surfmb = news ; Plane Surface(surfmb) = {llmbout,llmbin}; // Complete (otherwise is too complicated...)

//============================================================
//============================================================
If(Flag_Symmetry==0)
  // FULL MODEL
  // Rotation of Pi + duplication of all surfaces
  NN = #arcso[]-1 ;
  For N In {0:NN} // For simplicity (those lines would appear naturally when rotating Stator[0])
    arcso[]+= Rotate {{0, 0, 1}, {0, 0, 0}, Pi} { Duplicata{ Line{arcso[N]};} };
  EndFor

  Stator[]+= Rotate {{0, 0, 1}, {0, 0, 0}, Pi} { Duplicata{ Surface{Stator[0]};} };
  Rotor[]+= Rotate {{0, 0, 1}, {0, 0, 0}, Pi} { Duplicata{ Surface{Rotor[0]};} };
  Shaft[]+= Rotate {{0, 0, 1}, {0, 0, 0}, Pi} { Duplicata{ Surface{Shaft[0]};} };

  For N In {0:Ns/fac-1}
    Coiln[]+= Rotate {{0, 0, 1}, {0, 0, 0}, Pi} { Duplicata{ Surface{Coiln[N]};} };
    Coilp[]+= Rotate {{0, 0, 1}, {0, 0, 0}, Pi} { Duplicata{ Surface{Coilp[N]};} };
  EndFor

  surfmbstator[]+= Rotate {{0, 0, 1}, {0, 0, 0}, Pi} { Duplicata{ Surface{surfmbstator[0]};} };
  surfmbrotor[]+= Rotate {{0, 0, 1}, {0, 0, 0}, Pi} { Duplicata{ Surface{surfmbrotor[0]};} };

  NN = #AirgapRotorIn[]-1 ;
  For N In {0:NN}
    AirgapRotorIn[]+= Rotate {{0, 0, 1}, {0, 0, 0}, Pi} { Duplicata{ Surface{AirgapRotorIn[N]};} };
  EndFor
EndIf

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
// Physical regions
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

Physical Surface(STATOR) = {Stator[]} ;
Physical Surface(ROTOR)  = {Rotor[]} ;
Physical Surface(SHAFT)  = {Shaft[]} ;

NN = (!Flag_Symmetry)?Ns:Ns/2;
For N In {0:NN-1}
 Physical Surface(COILN+N) = {Coiln[N]} ;
 Physical Surface(COILP+N) = {Coilp[N]} ;
EndFor
Physical Surface(AIRGAP_STATOR) = {surfmbstator[]};
Physical Surface(AIRGAP_DUMMY)  = {surfmb[]}; // Not used
Physical Surface(AIRGAP_ROTOR)  = {surfmbrotor[]};
Physical Surface(AIRGAP_ROTORIN)= {AirgapRotorIn[]};

Physical Line(BND_STATOR) = arcso[] ;

If(!Flag_Symmetry)
  Physical Line(LMB_STATOR) = clmbs[] ;
  Physical Line(LMB_ROTOR)  = clmbr[] ;
EndIf
If(Flag_Symmetry)
  Physical Line(LMB_STATOR)   = clmbs[{0:#clmbs[]/2-1}] ;// half touching real geo
  Physical Line(LMB_ROTOR)    = clmbr[{0:#clmbr[]/2-1}] ;// half touching real geo
  Physical Line(LMB_ROTORAUX) = clmbr[{#clmbr[]/2:#clmbr[]-1}] ; //half for enaubling motion

  Physical Point(PMB_ROTORREF) = {pAxe} ;

  //Lines for symmetry link
  Physical Line(LSTATOR_RIGHT) = {ss1,clci1[0],lmbs[0]} ;
  Physical Line(LSTATOR_LEFT) = {ss2,clci1[(Ns/fac)%Ns],lmbs[1]} ;
  Physical Line(LROTOR_RIGHT) = {rr1,rr1_,cutshaft[1],lmbr[0]} ;
  Physical Line(LROTOR_LEFT) =  {rr2,rr2_,cutshaft[0],lmbr[1]};
EndIf

//-------------------------------------------------------------------------------
// For nice visualization
//-------------------------------------------------------------------------------
For N In {0:NN-1}
  linCoil[] += Boundary{Surface{Coiln[N], Coilp[N]};};
EndFor
linStator[] = CombinedBoundary{ Surface{Stator[]};};
linRotor[]  = CombinedBoundary{ Surface{Rotor[]}; };
Physical Line(LINE_NICEVIEW) = {linStator[],linRotor[],linCoil[]} ;

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

