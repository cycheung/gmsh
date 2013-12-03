Group{
  Omega = Region[7]; // Omega
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

FunctionSpace {
  { Name Hcurl_e; Type Form1;
    BasisFunction {
      // Ordre 1 Complet //
      { Name se   ; NameOfCoef ee   ; Function BF_Edge ; Support Region[Omega] ; Entity EdgesOf[All]; }
    }
  }
}

Formulation {
  { Name MaxwellA; Type FemEquation;
    Quantity {
      { Name e; Type Local;  NameOfSpace Hcurl_e; }
    }
    Equation {
      Galerkin { [ Dof{d e} , {d e} ];
        In Omega; Integration IOrder2; Jacobian JVol;  }
    }
  }

  { Name MaxwellB; Type FemEquation;
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
  { Name Maxwell_e ;
    System {
      { Name A ; NameOfFormulation MaxwellA ; }
      { Name B ; NameOfFormulation MaxwellB ; }
    }
    Operation {
      Generate[A] ; Print[A];
      Generate[B] ; Print[B];
    }
  }
}
