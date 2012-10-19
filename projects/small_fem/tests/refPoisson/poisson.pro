Group {
  Omega  = Region[ {9} ];
  Gamma1 = Region[ {5} ];
}

Jacobian {
  { Name Vol ; Case { { Region All ; Jacobian Vol ; } } }
}

Integration {
  { Name Int ;
    Case { {Type Gauss ;
            Case { { GeoElement Triangle    ; NumberOfPoints  4 ; }
                   { GeoElement Quadrangle  ; NumberOfPoints  4 ; }
                   { GeoElement Tetrahedron ; NumberOfPoints  4 ; }
                   { GeoElement Hexahedron  ; NumberOfPoints  6 ; }
                   { GeoElement Prism       ; NumberOfPoints  9 ; } }
           }
         }
  }
}

Constraint {
  { Name Dirich ;
    Case {
      { Region Gamma1 ; Value 0. ; }
    }
  }
}

FunctionSpace {
  { Name fs ; Type Form0 ;
    BasisFunction {
      { Name sn ; NameOfCoef vn ; Function BF_Node ;
        Support Omega ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf ; NameOfConstraint Dirich; }
    }
  }
}

Formulation {
  // Poisson equation with Dirichlet BC and source term
  { Name f ; Type FemEquation ;
    Quantity {
      { Name v1 ; Type Local ; NameOfSpace fs ; }
    }
    Equation {
      Galerkin { [ Dof{d v1} , {d v1} ] ; In Omega ; Jacobian Vol ; Integration Int ; }
      Galerkin { [ - 1. , {v1} ] ; In Omega ; Jacobian Vol ; Integration Int ; }
    }
  }
}

Resolution {
  { Name r ;
    System {
      { Name A1 ; NameOfFormulation f ; }
    }
    Operation {
      Generate[A1] ; Solve[A1] ; SaveSolution[A1] ;
    }
  }
}

PostProcessing {
  { Name p ; NameOfFormulation f  ;
    Quantity {
      { Name v1 ; Value { Term { [ {v1} ] ; In Omega ; Jacobian Vol; } } }
    }
  }
}

PostOperation{
  { Name p ; NameOfPostProcessing p;
    Operation {
      Print[ v1 , OnElementsOf Omega , File "poisson.pos" ] ;
    }
  }
}

