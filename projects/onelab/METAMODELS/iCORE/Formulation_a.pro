
/* -------------------------------------------------------------------------- 
   MagSta_a : Magnetic vector potential a formulation (2D) 
   -------------------------------------------------------------------------- */

Group {
  Dirichlet_a   = Region[ 1000 ] ;
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
  { Name a ;
    Case {
      { Region Dirichlet_a ; Value 0. ; }
    }
  }
}

FunctionSpace {

  { Name Hcurl_a ; Type Form1P ;
    BasisFunction {
      { Name se ; NameOfCoef ae ; Function BF_PerpendicularEdge ;
        Support Domain ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef ae ; EntityType NodesOf ; NameOfConstraint a ; }
    }
  }

}


Formulation {
  { Name MagSta ; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local ; NameOfSpace Hcurl_a ; }
    }
    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ] ; 
                 In Domain ; Jacobian JVol ; Integration I1 ; }

      Galerkin { [ hc[] , {d a} ] ; 
                 In Domain_M ; Jacobian JVol ; Integration I1 ; }

      Galerkin { [ -js[] , {a} ] ; 
                 In Domain_S ; Jacobian JVol ; Integration I1 ; }
    }
  }
}


Resolution {
  { Name MagSta ;
    System {
      { Name A ; NameOfFormulation MagSta ; }
    }
    Operation { 
      Generate[A] ; Solve[A] ; SaveSolution[A] ;
    }
  }
}


PostProcessing {
  { Name MagSta ; NameOfFormulation MagSta ;
    Quantity {
      { Name a ; Value { Local { [ CompZ[{a}] ]   ; In Domain ; Jacobian JVol ; } } }
      { Name b ; Value { Local { [ {d a} ]        ; In Domain ; Jacobian JVol ; } } }
      { Name a ; Value { Local { [ {a} ]          ; In Domain ; Jacobian JVol ; } } }
      { Name h ; Value { Local { [ nu[] * {d a} ] ; In Domain ; Jacobian JVol ; } 
                         Local { [ hc[] ]         ; In Domain_M ; Jacobian JVol ; } } }
    }
  }
}

eps = 1.e-5 ;
PostOperation {
  { Name MagSta ; NameOfPostProcessing MagSta;
    Operation {
      //Print[ a, OnElementsOf Domain, File "a.pos"] ;
      //Print[ a, OnElementsOf Domain, File > "a.pos", ChangeOfCoordinates {$X,-$Y,$Z}, ChangeOfValues {2*$Val0/2}] ;
      Print[ b, OnElementsOf Domain, File "b_a.pos" ] ;
      // Print[ h, OnElementsOf Domain, File "h_a.pos", Depth 0 ] ;
      //Print[ b, OnPlane {{-0.1,0,0}{0.1, 0, 0}{-0.1,0.1,0}} {60,30}, File "b_a_grid.pos" ] ;
      //Print[ b, OnLine {{-0.07,eps,0}{0.09, eps, 0}} {500}, File "b_a.cut" , Format Table ] ;
      //Print[ h, OnLine {{-0.06,eps,0}{-0.06,0.05,0}} {500}, File "h_a.cut" , Format Table ] ;
    }
  }
}
