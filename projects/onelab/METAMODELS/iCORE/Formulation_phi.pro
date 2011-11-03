/* --------------------------------------------------------------------------
   MagSta_phi : Magnetic scalar potential phi formulation 
   -------------------------------------------------------------------------- */

Group {
  Dirichlet_phi = Region[ {1000,1001} ] ;	
}

Integration {
  { Name I1 ;
    Case { 
      { Type Gauss ;
        Case { 
	  { GeoElement Triangle    ; NumberOfPoints  4 ; }
	  { GeoElement Quadrangle  ; NumberOfPoints  4 ; }
	}
      }
    }
  }
}

Constraint {
  { Name phi ;
    Case {
      { Region Dirichlet_phi ; Value 0. ; }
    }
  }
}

FunctionSpace {
  { Name Hgrad_phi ; Type Form0 ;
    BasisFunction {
      { Name sn ; NameOfCoef phin ; Function BF_Node ;
        Support Domain ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef phin ; EntityType NodesOf ; NameOfConstraint phi ; }
    }
  }
}

Formulation {
  { Name MagSta ; Type FemEquation ;
    Quantity { 
      { Name phi ; Type Local ; NameOfSpace Hgrad_phi ; }
    }
    Equation {
      Galerkin { [ - mu[] * Dof{d phi} , {d phi} ] ;  
                 In Domain ; Jacobian JVol ; Integration I1 ; }
      Galerkin { [ - mu[] * hc[] , {d phi} ] ;  
                 In Domain_M ; Jacobian JVol ; Integration I1 ; }
    }
  }
}

Resolution {
  { Name MagSta ;
    System {
      { Name A ; NameOfFormulation MagSta ; }
    }
    Operation { 
      Generate[A] ;
      //Print[A];
      Solve[A] ; SaveSolution[A] ; 
    }
  }
}

PostProcessing {
  { Name MagSta ; NameOfFormulation MagSta ;
    Quantity {
      { Name b   ; Value { Local { [ - mu[] * {d phi} ] ; In Domain ; Jacobian JVol ; } 
                           Local { [ - mu[] * hc[] ]    ; In Domain_M ; Jacobian JVol ; } } }
      { Name h   ; Value { Local { [ - {d phi} ]        ; In Domain ; Jacobian JVol ; } } }
      { Name phi ; Value { Local { [ {phi} ]            ; In Domain ; Jacobian JVol ; } } }
    }
  }
}

eps = 1.e-5 ;
PostOperation {
  { Name Magsta ; NameOfPostProcessing MagSta;
    Operation {
      //Print[ phi, OnElementsOf Domain, File "phi.pos"] ;
      Print[ b,   OnElementsOf Domain, File "b_phi.pos" ] ;
      //Print[ h,   OnElementsOf Domain, File "h_phi.pos", Depth 0 ] ;
      //Print[ b,   OnPlane {{-0.1,0,0}{0.1, 0, 0}{-0.1,0.1,0}} {60,30}, File "b_phi_grid.pos" ] ;
      //Print[ b,   OnLine {{-0.07,eps,0}{0.09, eps, 0}} {500}, File "b_phi.cut" , Format Table ] ;
      //Print[ h,   OnLine {{-0.06,eps,0}{-0.06,0.05,0}} {500}, File "h_phi.cut" , Format Table ] ;
      //Print[ phi, OnGrid {$A*Cos[$B]-0.05,$A*Sin[$B],0} { 0:0.05:0.001, 0:Pi:Pi/51, 0}, File "phi_polar.pos" ] ;
      //Echo [ "Plugin(Triangulate).iView = PostProcessing.NbViews-1;", File >>"phi_polar.pos" ];
      //Echo [ "Plugin(Triangulate).Run;", File >>"phi_polar.pos" ];
    }
  }
}