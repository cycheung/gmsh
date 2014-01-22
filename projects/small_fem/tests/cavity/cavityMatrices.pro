Group{
  Border = Region[5]; // Gamma
  Omega  = Region[7]; // Omega
}

Jacobian {
  { Name JVol ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }
}

Integration {
  { Name IOrder2 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Line ; NumberOfPoints  2 ; }
          { GeoElement Triangle ; NumberOfPoints 3 ; }
          { GeoElement Tetrahedron ; NumberOfPoints 4 ; }
        }
      }
    }
  }

  { Name IOrder4 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Line ; NumberOfPoints  3 ; }
          { GeoElement Triangle ; NumberOfPoints 6 ; }
          { GeoElement Tetrahedron ; NumberOfPoints 15 ; }
        }
      }
    }
  }
}

Constraint {
  { Name Dirichlet_e;
    Case {
      { Region Border ; Value 0 ; }
    }
  }
}

FunctionSpace {
  { Name Hcurl_e ; Type Form1;
    BasisFunction {
      { Name se ; NameOfCoef ee ; Function BF_Edge ; Support Region[{Omega,Border}] ; Entity EdgesOf[All]; }
    }

    Constraint {
      { NameOfCoef ee   ; EntityType EdgesOf  ; NameOfConstraint Dirichlet_e; }
    }
  }
}

Formulation {
  { Name Maxwell_A; Type FemEquation;
    Quantity {
      { Name e; Type Local;  NameOfSpace Hcurl_e; }
    }
    Equation {
      Galerkin { [ Dof{d e} , {d e} ];
        In Omega; Integration IOrder2; Jacobian JVol;  }
    }
  }

  { Name Maxwell_B; Type FemEquation;
    Quantity {
      { Name e; Type Local;  NameOfSpace Hcurl_e; }
    }
    Equation {
      Galerkin { [ Dof{e} , {e} ];
        In Omega; Integration IOrder4; Jacobian JVol;  }
    }
  }
}

Resolution {
  { Name Maxwell_A ;
    System {
      { Name A ; NameOfFormulation Maxwell_A ; }
    }
    Operation {
      Generate[A] ; Print[A] ;
    }
  }

  { Name Maxwell_B ;
    System {
      { Name B ; NameOfFormulation Maxwell_B ; }
    }
    Operation {
      Generate[B] ; Print[B] ;
    }
  }
}
