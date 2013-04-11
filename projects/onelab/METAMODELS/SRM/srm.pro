
Include "srm_data.pro";

/*
These variables are pointles when the model is called from a python metamodel

DefineConstant[ ComputeCommand = {"-solve -v 1", Path "getdp/9"} ];
// with "-solve -v 1 -v2" a mesh-based Gmsh output is issued, i.e.
// all view are appended to the same file and appear superimposed in the GUI 
DefineConstant[ ResolutionChoices = {"MagSta", Path "getdp/1"} ];
DefineConstant[ PostOperationChoices = {"Torque", Path "getdp/2"} ];
*/

/*
Re-declaration of onelab parameter, to make the .pro file selfconsistent
The parameters are already declared in the main python script
*/
DefineConstant[ POSITION = {30., Path "IO/1", Label "Rotor Position [deg]"}];
DefineConstant[ MST = {0., Path "IO/2", Label "Torque (MST) [Nm]"}];
DefineConstant[ VWP = {0., Path "IO/3", Label "Torque (VWP) [Nm]"}];
DefineConstant[ IR = {  1, Path "IO/4", Label "Current phase R [A]"} ];
DefineConstant[ IS = {  0, Path "IO/5", Label "Current phase S [A]"} ];
DefineConstant[ IT = {  0, Path "IO/6", Label "Current phase T [A]"} ];

angleR = POSITION ;// angular rotor position


Flag_NL  = 0  ; // 1 = BH curve for Iron (nonlinear case)
Flag_Cir = 0  ; // 1 = Circuit for imposing source

Group {
  AirgapStator = #AIRGAP_STATOR;
  AirgapRotor  = #{AIRGAP_ROTOR};
  AirgapRotorIn  = #{AIRGAP_ROTORIN};
  Shaft  = #SHAFT;
  Stator = #STATOR; Rotor  = #ROTOR;

  Outer =  #BND_STATOR ;
  LinStatorRight = #LSTATOR_RIGHT ;
  LinStatorLeft  = #LSTATOR_LEFT;
  LinRotorRight  = #LROTOR_RIGHT;
  LinRotorLeft   = #LROTOR_LEFT;

  // Moving band
  LinMBStator  = #LMB_STATOR;
  LinMBRotor_1 = #LMB_ROTOR;
  LinMBRotor_2 = #LMB_ROTORAUX;
  LinMBRotorAux = Region[{LinMBRotor_2}];
  LinMBRotor   = Region[{LinMBRotor_1,LinMBRotor_2}];
  PntMBRotorCen = #PMB_ROTORREF;

  FactorSym = (!Flag_Symmetry)? 1.:2.;
  MB = MovingBand2D[9999, LinMBStator, LinMBRotor, FactorSym];
  RotorMoving = #{Rotor, AirgapRotor, AirgapRotorIn, Shaft, LinMBRotorAux}; // Use in Change of Coordinates

  Air = #{AirgapStator, AirgapRotor, AirgapRotorIn, Shaft, MB};

  For i In {1:Ns/FactorSym}
    IndP~{i} = Region[{(COILP+i-1)}] ;
    IndN~{i} = Region[{(COILN+i-1)}] ;
    Inds  += Region[{IndP~{i},IndN~{i}}] ;
  EndFor

  //Phase_1 = #{IndP_1,IndN_1,IndP_4,IndN_4};
  //Phase_2 = #{IndP_2,IndN_2,IndP_5,IndN_5};
  //Phase_3 = #{IndP_3,IndN_3,IndP_6,IndN_6};
  For ip In {1:Ns/2} //3 phases
    For i In {1:2/FactorSym}
      Phase~{ip} += Region[{IndP~{(i-1)*3+ip}, IndN~{(i-1)*3+ip}}];
    EndFor
  EndFor

  For ip In {1:Ns/2}
    IndIdirP += Region[{IndP~{ip}}];
    IndIdirN += Region[{IndN~{ip}}];
    For i In {1:2/FactorSym-1}
      IndIdirP += Region[{IndN~{i*3+ip}}];
      IndIdirN += Region[{IndP~{i*3+ip}}];
    EndFor
  EndFor

  DomainS = Region[{Phase_1, Phase_2, Phase_3}];
  DomainB = Region[{}];
  If(Flag_Cir)
    DomainS = Region[{}];
    DomainB = Region[{Phase_1, Phase_2, Phase_3}];
  EndIf

  DomainL  = Region[{Air, Inds, Stator,Rotor}];
  DomainNL = Region[{}];
  If(Flag_NL)
    DomainL  = Region[{Stator,Rotor}];
    DomainNL = Region[{Stator,Rotor}];
  EndIf

  Domain = Region[{ DomainL, DomainNL }];

  If(Flag_Cir)
    // Dummy numbers for circuit definition
    R_1 = #551 ;
    R_2 = #552 ;
    R_3 = #553 ;

    E_1 = #771 ;
    E_2 = #772 ;
    E_3 = #773 ;

    //Resistance_Cir = #{R_1,R_2,R_3};
    //Source_Cir     = #{E_1,E_2,E_3};
    Resistance_Cir = #{R_1};
    Source_Cir     = #{E_1};
    DomainZt_Cir   = #{Resistance_Cir, Source_Cir};
  EndIf

  DomainKin = #1234 ; // Dummy region number for mechanical equation

}

Function {
  mu0 = 4.e-7 * Pi ;
  murIron = 1000 ;
  sigmaIron = 5e6 ;

  nu [ #{Air, Inds} ] = 1. / mu0 ;
  nuIron = 1./(mu0*murIron) ; // For the linear case

  nu[ #{Stator,Rotor} ] = nuIron ;
  sigma[ #{Stator,Rotor} ] = sigmaIron ;

  // h[]-vector and dhdb[]-tensor as function of b-vector for nonlinear materials
  // Analytical BH-curve
  nu_1a[] = 100. + 10. * Exp[1.8*$1] ;
  dnudb2_1a[] = 18. * Exp[1.8*$1] ;
  h_1a[] = nu_1a[SquNorm[$1]] * $1 ;
  dhdb_1a[] = TensorDiag[1,1,1] * nu_1a[SquNorm[$1]#1] + 2*dnudb2_1a[#1] * SquDyadicProduct[$1]  ;

  h [ DomainNL ]    = h_1a[$1];
  dhdb [ DomainNL ] = dhdb_1a[$1];

  /*
  // linear law => testing purpose
  h [ Domain_NonLin ] = nuIron * $1 ;
  dhdb [ Domain_NonLin ] = nuIron * TensorDiag[1,1,1] ;
  */

  If(Flag_Cir)
    Vdc = 220 * Sqrt[2]/FactorSym ;
    Resistance[#DomainB]  = 1 ;
    Resistance[#{Resistance_Cir}]  = 1 ;
  EndIf

  // Parameters for Nonlinear Loop
  NL_Eps=1e-3; NL_Relax=1; NL_NbrMax=100;

  // fixed rotor position or inital position (in rad) in case of rotation
  th0 = angleR*Pi/180 ;
  // end angle (in rad)
  th1 = (angleR+180)*Pi/180 ; 

  // imposed speed and current profile per phase
  rpm = 5000 ; // in rpm
  Om = rpm/60*2*Pi ; // in rad/s

  nth = 10 ; // 100 number of angle and time steps
  dth = (th1-th0)/nth ; // angle step (in rad)
  dt = dth/Om ; // time step

  dt  = 1e-3/2 ; //nth = 200 ;
  t0 = 0 ;  tmax = nth * dt ;

  Nw = 226 ; // Number of turns per coil
  Idir[#{IndIdirP}] =  Nw ;
  Idir[#{IndIdirN}] = -Nw ;
  NbrWires_Area[] = Nw/SurfaceArea[]{COILN} ;

  // supply at fixed position
  Freq = 50 ;
  Omega = 2*Pi*Freq ;
  T = 1/Freq ;

  If(!Flag_Cir)
    Imax = 1 ;
    // UI[] = F_Cos_wt_p[]{2*Pi*Freq, 0.};
    // I_1[Phase_1] = UI[] ;
    // I_2[Phase_2] = 0. ;
    // I_3[Phase_3] = 0. ;

    I_1[Phase_1] = IR ;
    I_2[Phase_2] = IS ;
    I_3[Phase_3] = IT ;
    Iphase[#Phase_1] = Imax * Idir[] * I_1[] ;
    Iphase[#Phase_2] = Imax * Idir[] * I_2[] ;
    Iphase[#Phase_3] = Imax * Idir[] * I_3[] ;
    is[] = Iphase[];
  EndIf

  // Change of coordinates
  RotatePZ[] = Rotate[ Vector[$X,$Y,$Z], 0, 0, $1 ] ;
  // Maxwell stress tensor
  T_max[] = ( SquDyadicProduct[$1] - SquNorm[$1] * TensorDiag[0.5,0.5,0.5] ) / mu0 ;

  // Kinematics
  Inertia = 8.3e-3 ;
  Friction[] = 0 ;

  //Fmag[] =  1000*F_Cos_wt_p[]{2*Pi*Freq, 0.}; // testing
  Fmag[] = #55 ; // Computed in postprocessing (See "Store 55" below)
}



Constraint {
  { Name MVP_2D ; // Magnetic Vector Potential - a
    Case {
      { Region Outer ; Value 0. ; }

      If(Flag_Symmetry)
        { Region LinStatorLeft; SubRegion Outer; Type Link;
          RegionRef LinStatorRight; SubRegionRef Outer; Coefficient -1.; Function RotatePZ[Pi]; }
        { Region LinRotorLeft; SubRegion PntMBRotorCen; Type Link;
          RegionRef LinRotorRight; SubRegionRef PntMBRotorCen;  Coefficient -1.; Function RotatePZ[Pi]; }
        { Region LinMBRotor_2; SubRegion LinMBRotor_1; Type Link;
          RegionRef LinMBRotor_1; SubRegionRef LinMBRotor_2;  Coefficient -1.; Function RotatePZ[Pi]; }
      EndIf
    }
  }

  If(Flag_Cir)
  If(!Flag_Symmetry)
  { Name ElectricalCircuit ; Type Network ;
    Case Circuit_1 {
      { Region E_1 ;    Branch {100,101} ; }
      { Region R_1 ;    Branch {101,102} ; }
      { Region IndP_1 ; Branch {102,103} ; }
      { Region IndP_4 ; Branch {104,103} ; }
      { Region IndN_1 ; Branch {105,104} ; }
      { Region IndN_4 ; Branch {105,100} ; }
    }
    /*
    Case Circuit_2 {
      { Region E_2 ;    Branch {200,201} ; }
      { Region R_2 ;    Branch {201,202} ; }
      { Region IndP_2 ; Branch {202,203} ; }
      { Region IndP_5 ; Branch {204,203} ; }
      { Region IndN_2 ; Branch {205,204} ; }
      { Region IndN_5 ; Branch {205,200} ; }
    }
    Case Circuit_3 {
      { Region E_3 ;    Branch {300,301} ; }
      { Region R_3 ;    Branch {301,302} ; }
      { Region IndP_3 ; Branch {302,303} ; }
      { Region IndP_6 ; Branch {304,303} ; }
      { Region IndN_3 ; Branch {305,304} ; }
      { Region IndN_6 ; Branch {305,300} ; }
    }
    */
  }
  EndIf
  If(Flag_Symmetry)
  { Name ElectricalCircuit ; Type Network ;
    Case Circuit_1 {
      { Region E_1 ;    Branch {100,101} ; }
      { Region R_1 ;    Branch {101,102} ; }
      { Region IndP_1 ; Branch {102,103} ; }
      { Region IndN_1 ; Branch {100,103} ; }
    }
    /*
    Case Circuit_2 {
      { Region E_2 ;    Branch {200,201} ; }
      { Region R_2 ;    Branch {201,202} ; }
      { Region IndP_2 ; Branch {202,203} ; }
      { Region IndN_2 ; Branch {200,203} ; }
    }
    Case Circuit_3 {
      { Region E_3 ;    Branch {300,301} ; }
      { Region R_3 ;    Branch {301,302} ; }
      { Region IndP_3 ; Branch {302,303} ; }
      { Region IndN_3 ; Branch {300,303} ; }
    }
    */
  }
  EndIf

  { Name Voltage_Cir ;
    Case {
      { Region E_1 ; Value Vdc ;}
    }
  }
  { Name Current_Cir ;
    Case {
    }
  }

  // Stranded conductors in Domain B
  { Name Voltage_2D ;
    Case {
    }
  }

  { Name Current_2D ;
    Case {
      { Region Phase_2 ; Value 0. ;}
      { Region Phase_3 ; Value 0. ;}
    }
  }
 EndIf

 //Kinetics
 { Name CurrentPosition ;
   Case {
     { Region DomainKin ; Type Init ; Value 0.#66 ; }  
   }
 }

 { Name CurrentVelocity ;
    Case {
      { Region DomainKin ; Type Init ; Value 0. ; }
    }
  }

}


Dir = "res/";
ExtGmsh     = Str[Sprintf("%g_dt%g.pos",FactorSym, dt) ];
ExtGnuplot  = Str[Sprintf("%g_dt%g.dat",FactorSym, dt) ];

cleanT = StrCat["rm -rf ", StrCat[Dir, StrCat["T", ExtGnuplot]] ];
cleanP = StrCat["rm -rf ", StrCat[Dir, StrCat["P", ExtGnuplot]] ];
cleanV = StrCat["rm -rf ", StrCat[Dir, StrCat["V", ExtGnuplot]] ];

Include "MagSta2D_MB.pro"

PostOperation ab_lastTS UsingPost MagSta_a_2D {
  Print[ a, OnElementsOf #{Domain, LINE_NICEVIEW}, File StrCat[Dir, StrCat["a",ExtGmsh]], AppendTimeStepToFileName] ;
  //Print[ b, OnElementsOf #{Domain, LINE_NICEVIEW}, File StrCat[Dir, StrCat["b",ExtGmsh]], AppendTimeStepToFileName] ; //LastTimeStepOnly
}

PostOperation FieldLines UsingPost MagSta_a_2D {
  Print[ a, OnElementsOf #{Domain, LINE_NICEVIEW}, File "a.pos"] ;
}

PostOperation Torque UsingPost MagSta_a_2D {
  Print[ T_Maxwell[AirgapStator], OnGlobal, Format TimeTable, File > StrCat[Dir, StrCat["T", ExtGnuplot]],
         LastTimeStepOnly, Store 55, SendToServer "IO/2MST" ] ;
  Print[ T_vw, OnRegion NodesOf[LinMBStator], Format RegionValue, File StrCat[Dir, StrCat["Tvw", ExtGnuplot]], LastTimeStepOnly, SendToServer "IO/3VWP" ] ;
  //Print[ T_Maxwell[AirgapRotor], OnGlobal, Format TimeTable, File StrCat[Dir, StrCat["Tr", ExtGnuplot]], LastTimeStepOnly ] ;
  //Print[ T_vwR, OnRegion NodesOf[LinMBRotor], Format RegionValue, File StrCat[Dir, StrCat["Tvwr", ExtGnuplot]], LastTimeStepOnly ] ;
}

PostOperation flux_UI UsingPost MagSta_a_2D {
  Print[ a_av[#Phase_1], OnGlobal , Format TimeTable , File StrCat["flux1", ExtGnuplot]] ;
  Print[ a_av[#Phase_2], OnGlobal , Format TimeTable , File StrCat["flux2", ExtGnuplot]] ;
  Print[ a_av[#Phase_3], OnGlobal , Format TimeTable , File StrCat["flux3", ExtGnuplot]] ;

  If(Flag_Cir)
    Print[ U, OnRegion Source_Cir, Format Table, File StrCat[Dir, StrCat["Uc",ExtGnuplot]] ] ;
    Print[ I, OnRegion Resistance_Cir, Format Table, File StrCat[Dir, StrCat["Ic",ExtGnuplot]] ] ;
  EndIf
}

PostOperation MapMec UsingPost Mechanical { // angular rotor speed (rad/s) and position (rad)

  Print[ P, OnRegion DomainKin, File > StrCat[Dir, StrCat["P", ExtGnuplot]], Format Table, Store 77,
         LastTimeStepOnly, SendToServer "IO/1Position"] ;
  Print[ V, OnRegion DomainKin, File > StrCat[Dir, StrCat["V", ExtGnuplot]], Format Table,
         LastTimeStepOnly, SendToServer "Parameters/4Velocity"] ;
}
