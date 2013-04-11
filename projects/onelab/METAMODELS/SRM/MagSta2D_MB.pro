Group{
  DefineGroup[ DomainZt_Cir ] ;
  DefineGroup[ DomainKin ] ;
}

Function{
  DefineFunction[ Friction, Fmag ] ;
  DefineVariable[ Inertia ];
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

Jacobian {
  { Name J ; Case { { Region All ; Jacobian Vol ; } } }
}

Integration {
  { Name I ; Case { {
        Type Gauss ;
        Case {
	  { GeoElement Triangle    ; NumberOfPoints  6 ; }
	  { GeoElement Quadrangle  ; NumberOfPoints  4 ; }
	  { GeoElement Line        ; NumberOfPoints  13 ; }
	} } }
  }
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
FunctionSpace {

  { Name Hcurl_a_2D ; Type Form1P ;
    BasisFunction {
      { Name se ; NameOfCoef ae ; Function BF_PerpendicularEdge ;
        Support #{Domain,LinMBRotorAux} ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef ae  ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
    }
  }

  { Name Hregion_i_2D ; Type Vector ; // Stranded conductors
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_RegionZ ;
        Support DomainB ; Entity DomainB ; }
    }
    GlobalQuantity {
      { Name Ib ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Ub ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Ub ; EntityType Region ; /*NameOfConstraint Voltage_2D ;*/ }
      { NameOfCoef Ib ; EntityType Region ; /*NameOfConstraint Current_2D ;*/ }
    }
  }

  { Name Hregion_Z ; Type Scalar ; // Circuit equations
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_Region ;
        Support DomainZt_Cir ; Entity DomainZt_Cir ; }
    }
    GlobalQuantity {
      { Name Iz ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Uz ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Uz ; EntityType Region ; /*NameOfConstraint Voltage_Cir ;*/ }
      { NameOfCoef Iz ; EntityType Region ; /*NameOfConstraint Current_Cir ;*/ }
    }
  }

  { Name Position ; Type Scalar ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_Region ;
        Support DomainKin ; Entity DomainKin ; }
    }
    GlobalQuantity {
      { Name P ; Type AliasOf  ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef P ; EntityType Region ; NameOfConstraint CurrentPosition ; }
    }
  }

  { Name Velocity ; Type Scalar ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_Region ;
        Support DomainKin ; Entity DomainKin ; } }
    GlobalQuantity {
      { Name V ; Type AliasOf  ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef V ; EntityType Region ; NameOfConstraint CurrentVelocity ; }
    }
  }

}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

Formulation {
  //--------------------------------------------------------------------------
  // Magnetostatics
  //--------------------------------------------------------------------------
  { Name MagSta_a_2D ; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local ; NameOfSpace Hcurl_a_2D ; }

      { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
      { Name Ub ; Type Global ; NameOfSpace Hregion_i_2D [Ub] ; }
      { Name Ib ; Type Global ; NameOfSpace Hregion_i_2D [Ib] ; }

      { Name Uz ; Type Global ; NameOfSpace Hregion_Z [Uz] ; }
      { Name Iz ; Type Global ; NameOfSpace Hregion_Z [Iz] ; }
    }

    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ] ;
        In DomainL ; Jacobian J ; Integration I ; } // If(Flag_NL) Domain = DomainL ;

      Galerkin { [ h[{d a}] , {d a} ] ;
        In DomainNL ; Jacobian J ; Integration I ; } // If(!Flag_NL) DomainNL = #{};
      Galerkin { JacNL [ dhdb[{d a}] * Dof{d a} , {d a} ] ;
        In DomainNL ; Jacobian J ; Integration I ; }
      Galerkin {  [  0*Dof{d a} , {d a} ]  ; // DO NOT REMOVE!!! - Keeping track of Dofs in auxiliary line of MB if Symmetry=1
        In LinMBRotorAux; Jacobian J; Integration I; }

      If(!Flag_Cir)
        Galerkin { [ -is[]/SurfaceArea[] * Vector[0, 0, 1]  , {a} ] ;
         In DomainS ; Jacobian J ; Integration I ; } //If (!Flag_Cir) DomainS = Inds
      EndIf

      If(Flag_Cir)
        Galerkin { [ -NbrWires_Area[] * Dof{ir} , {a} ] ;
          In DomainB ; Jacobian J ; Integration I ; } //If (!Flag_Cir) DomainB = Inds
        Galerkin { DtDof [ AxialLength * NbrWires_Area[]* Dof{a} , {ir} ] ;
          In DomainB ; Jacobian J ; Integration I ; }

        GlobalTerm { [ Dof{Ub}  , {Ib} ] ; In DomainB ; }
        GlobalTerm { [ Resistance[]  * Dof{Ib} , {Ib} ] ; In DomainB ; }
        GlobalTerm { NeverDt[ Dof{Uz}        , {Iz} ] ; In Resistance_Cir ; }
        GlobalTerm { NeverDt[ Resistance[{Uz}] * Dof{Iz} , {Iz} ] ; In Resistance_Cir ; }

        GlobalTerm {  [ 0*Dof{Uz}        , {Iz} ] ; In DomainZt_Cir ; }

        GlobalEquation {
          Type Network ; NameOfConstraint ElectricalCircuit ;
          { Node {Ib}; Loop {Ub}; Equation {Ub}; In DomainB ; }
          { Node {Iz}; Loop {Uz}; Equation {Uz}; In DomainZt_Cir ; }
         }
      EndIf
    }
  }

  //--------------------------------------------------------------------------
  // Mechanics
  //--------------------------------------------------------------------------
  { Name Mechanical ; Type FemEquation ;
    Quantity {
      { Name V ; Type Global ; NameOfSpace Velocity [V] ; } // velocity
      { Name P ; Type Global ; NameOfSpace Position [P] ; } // position
    }
    Equation {
      GlobalTerm { DtDof [ Inertia * Dof{V} , {V} ] ; In DomainKin ; }
      GlobalTerm { [ Friction[] * Dof{V} , {V} ] ; In DomainKin ; }
      GlobalTerm { [             -Fmag[] , {V} ] ; In DomainKin ; }

      GlobalTerm { DtDof [ Dof{P} , {P} ] ; In DomainKin ; }
      GlobalTerm {       [-Dof{V} , {P} ] ; In DomainKin ; }
    }
  }

}


//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
Resolution {

  { Name MagSta ;
    System {
      { Name A ; NameOfFormulation MagSta_a_2D ; }
    }
    Operation {
      ChangeOfCoordinates [NodesOf[RotorMoving], RotatePZ[th0]] ;
      InitMovingBand2D[MB] ; MeshMovingBand2D[MB] ;
      //SaveMesh[A, Region[{Domain}], "Domain.msh"]; // Checking mesh (e.g. moving band)
      Generate[A] ; Solve[A] ; SaveSolution[A] ;
      PostOperation[FieldLines] ;
      PostOperation [Torque] ;
    }
  }

  { Name  MagSta_time ; // Imposed movement
    System {
      { Name A ; NameOfFormulation MagSta_a_2D ; }
    }
    Operation {
      SystemCommand[cleanT];
      ChangeOfCoordinates [NodesOf[RotorMoving], RotatePZ[th0]] ;
      InitMovingBand2D[MB] ; MeshMovingBand2D[MB] ;
      InitSolution[A];

      TimeLoopTheta [t0, tmax, dt, 1.] {
        Generate[A] ; Solve[A] ; SaveSolution[A] ;
        PostOperation [ab_lastTS] ;
        PostOperation [Torque] ;
        ChangeOfCoordinates [NodesOf[RotorMoving], RotatePZ[dth]] ;
        MeshMovingBand2D[MB] ;
      }
    }
  }

  { Name  Mechanics ;
    System {
      { Name M ; NameOfFormulation Mechanical ; }
    }
    Operation {
      SystemCommand[cleanP]; SystemCommand[cleanV];
      InitSolution[M] ; SaveSolution[M] ;
      TimeLoopTheta[0.,1000*1e-3,1e-3, 1.]{
        Generate[M] ; Solve[M] ; SaveSolution[M] ;
        PostOperation[MapMec] ;
      }
    }
  }

  { Name  MagSta_Kin ;
    System {
      { Name A ; NameOfFormulation MagSta_a_2D ; }
      { Name M ; NameOfFormulation Mechanical ; }
    }
    Operation {
      SystemCommand[cleanT]; SystemCommand[cleanP]; SystemCommand[cleanV];

      ChangeOfCoordinates [ NodesOf[RotorMoving], RotatePZ[th0] ] ; // Initial position (supposing initial mesh with angleR=0)
      InitMovingBand2D[MB] ; MeshMovingBand2D[MB] ;

      InitSolution[A] ; SaveSolution[A] ;
      InitSolution[M] ; SaveSolution[M] ;

      TimeLoopTheta[t0, tmax, dt, 1.]{
	Generate[A] ; Solve[A] ;  SaveSolution[A] ;
        PostOperation[ab_lastTS] ;
        PostOperation[Torque] ;

        Generate[M] ; Solve[M] ; SaveSolution[M] ;
        PostOperation[MapMec] ;

        ChangeOfCoordinates [ NodesOf[RotorMoving], RotatePZ[#77-#66] ] ;
        Evaluate[ #77#66 ] ; //Keep track of previous angular position
        MeshMovingBand2D[MB] ;
      }
    }
  }

  { Name MagSta_NonLin ;
    System {
      { Name A ; NameOfFormulation MagSta_a_2D ; }
    }
    Operation {
      ChangeOfCoordinates [NodesOf[RotorMoving], RotatePZ[th0]] ;
      InitMovingBand2D[MB] ; MeshMovingBand2D[MB] ;
      InitSolution[A] ;
      IterativeLoop[NL_NbrMax, NL_Eps, NL_Relax] {
        GenerateJac[A]; SolveJac[A];
      }
      SaveSolution[A] ;
      PostOperation [ab_lastTS] ;
      PostOperation [Torque] ;
    }
  }

  { Name MagSta_NonLin_time ;
    System {
      { Name A ; NameOfFormulation MagSta_a_2D ; }
    }
    Operation {
      ChangeOfCoordinates [NodesOf[RotorMoving], RotatePZ[th0]] ;
      InitMovingBand2D[MB] ; MeshMovingBand2D[MB] ;
      InitSolution A;
      TimeLoopTheta [t0, tmax, dt, 1.] {
        IterativeLoop[NL_NbrMax, NL_Eps, NL_Relax] {
          GenerateJac A; SolveJac A;
        }
        SaveSolution A ;
        PostOperation [ab_lastTS] ;
        ChangeOfCoordinates [NodesOf[RotorMoving], RotatePZ[dth]] ;
        MeshMovingBand2D[MB] ;
      }
    }
  }

}


PostProcessing {

  { Name MagSta_a_2D ; NameOfFormulation MagSta_a_2D ;
    Quantity {
      { Name a ; Value { Local { [ CompZ[{a}] ]   ; In Domain ; Jacobian J ; } } }
      { Name b ; Value { Local { [ {d a} ]        ; In Domain ; Jacobian J ; } } }
      { Name h ; Value { Local { [ nu[] * {d a} ] ; In DomainL ; Jacobian J ; }
                         Local { [ h[{d a}] ] 	  ; In DomainNL ; Jacobian J; } } }
      { Name a_av ; Value {
          Integral { [ Idir[]*CompZ[{a}]/SurfaceArea[]* AxialLength ] ;
            In DomainS ; Jacobian J; Integration I; } } }

      { Name F_Maxwell  ; Value { Integral { [ AxialLength * T_max[{d a}] * Normal[] ] ;
           In DomainL ; Jacobian J ; Integration I ; } } }
      { Name T_Maxwell ; Value { // Torque from Maxwell stress tensor
          Integral { [ CompZ [ XYZ[] *^ (T_max[{d a}] * XYZ[] ) ] * 2 * Pi * AxialLength / SurfaceArea[] ] ;
            In DomainL ; Jacobian J ; Integration I; } } }

      { Name F_vw ; Value { // Magnetic force via virtual works
	  Integral { Type Global ;
            [ 0.5*nu[]*VirtualWork [{d a}] * AxialLength ];
		In DomainL; Jacobian J ; Integration I ; } } }

      { Name T_vw ; Value { // Torque from virtual work force
	Integral { Type Global ;
		[ CompZ[ 0.5 * nu[] * XYZ[] /\ VirtualWork[{d a}] ] * AxialLength ];
		In ElementsOf[AirgapStator, OnOneSideOf LinMBStator] ; Jacobian J ; Integration I; } } }
      { Name T_vwR ; Value { // Torque from virtual work force
	Integral { Type Global ;
		[ CompZ[-0.5 * nu[] * XYZ[] /\ VirtualWork[{d a}] ] * AxialLength ];
                In ElementsOf[AirgapRotor, OnOneSideOf LinMBRotor] ; Jacobian J ; Integration I; } } }
      If(Flag_Cir)
        { Name U ; Value { Term { [ {Ub} ]  ; In DomainB ; }
                           Term { [ {Uz} ]  ; In DomainZt_Cir ; } } }
        { Name I ; Value { Term { [ {Ib} ]  ; In DomainB ; }
                           Term { [ {Iz} ]  ; In DomainZt_Cir ; } } }
      EndIf
    }
  }

 { Name Mechanical ; NameOfFormulation Mechanical ;
   PostQuantity {
     { Name P  ; Value { Term { [ {P} ]  ; In DomainKin ; } } } //Position
     { Name V  ; Value { Term { [ {V} ]  ; In DomainKin ; } } } //Velocity
     { Name Vrpm  ; Value { Term { [ {V}*30/Pi ]  ; In DomainKin ; } } } //Velocity in rpm
   }
 }

}
