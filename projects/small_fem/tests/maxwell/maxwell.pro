Group{
  GammaS = Region[5]; // Source
  GammaC = Region[6]; // Conductor -- Reflection
  Omega  = Region[7]; // Omega
}


Function{
  k = 1;
}


Jacobian {
  { Name JVol ;
    Case { 
      { Region All ; Jacobian Vol ; }
    }
  }
}


Integration {
  { Name I1 ;
    Case { 
      { Type Gauss ;
        Case { 
          { GeoElement Point ; NumberOfPoints  1 ; }
          { GeoElement Line ; NumberOfPoints  4 ; }
          { GeoElement Triangle ; NumberOfPoints 12 ; }
          { GeoElement Quadrangle ; NumberOfPoints 7 ; }
          { GeoElement Tetrahedron ; NumberOfPoints 15 ; }
          { GeoElement Hexahedron ; NumberOfPoints 34 ; }
        }
      }
    }
  }
}


Constraint{
  { Name Dirichlet_e ;
    Case {      
      { Region GammaC ; Value  0. ; }
      { Region GammaS ; Value -1. ; }
    }
  }
}


FunctionSpace {
  { Name Hcurl_e; Type Form1;
    BasisFunction {
      { Name se; NameOfCoef ee; Function BF_Edge; Support Region[{Omega,GammaS,GammaC}] ; Entity EdgesOf[All]; }
    }
    Constraint {
      { NameOfCoef ee; EntityType EdgesOf ; NameOfConstraint Dirichlet_e; }
    }
  }
}


Formulation {
  { Name Maxwell_e; Type FemEquation; 
    Quantity { 
      { Name e; Type Local;  NameOfSpace Hcurl_e; }
    }
    Equation {
      Galerkin { [ Dof{d e} , {d e} ]; 
                 In Omega; Integration I1; Jacobian JVol;  }

      Galerkin { [ -k^2 * Dof{e} , {e} ]; 
                 In Omega; Integration I1; Jacobian JVol;  }      
    }
  }
}


Resolution {
  { Name Maxwell_e ;
    System {
      { Name A ; NameOfFormulation Maxwell_e ; } //Type Complex; }
    }
    Operation { 
      Generate[A] ; Solve[A] ; SaveSolution[A] ;
    }
  }
}


PostProcessing {
  { Name Maxwell_e ; NameOfFormulation Maxwell_e ;
    Quantity {
      { Name e ; 
	Value { Local { [ {e} ] ; In Omega; Jacobian JVol ; } } }
    }
  }
}


PostOperation {
  { Name Maxwell_e ; NameOfPostProcessing Maxwell_e; 
    Operation {
      Print[ e, OnElementsOf Omega, File "e.pos"] ;
    }
  }
}
