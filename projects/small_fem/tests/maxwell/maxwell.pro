Group{
  GammaS = Region[5]; // Source
  GammaC = Region[6]; // Conductor -- Reflection
  Omega  = Region[7]; // Omega
}


Function{
  k = 5;
}


Jacobian {
  { Name JVol ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }

  { Name JSur ;
    Case {
      { Region All ; Jacobian Sur ; }
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
      { Region GammaC ; Value Vector [0., 0., 0.] ; }
      { Region GammaS ; Value Vector [0., 0., 1.] ; }
    }
  }
}


FunctionSpace {
  { Name Hcurl_e; Type Form1;
    BasisFunction {
      // Ordre 1 Complet //
      { Name se   ; NameOfCoef ee   ; Function BF_Edge_1E ; Support Region[{Omega,GammaS,GammaC}] ; Entity EdgesOf[All]; }
      { Name se2e ; NameOfCoef we2e ; Function BF_Edge_2E ; Support Region[{Omega,GammaS,GammaC}] ; Entity EdgesOf[All]; }
    }
  }

  { Name HcurlLS_e; Type Form1;
    BasisFunction {
      // Ordre 1 Complet //
      { Name se   ; NameOfCoef ee  ; Function BF_Edge_1E; Support Region[{GammaS}] ; Entity EdgesOf[All]; }
      { Name se2e ; NameOfCoef we2e; Function BF_Edge_2E; Support Region[{GammaS}] ; Entity EdgesOf[All]; }
    }
  }

  { Name HcurlLC_e; Type Form1;
    BasisFunction {
      // Ordre 1 Complet //
      { Name se   ; NameOfCoef ee  ; Function BF_Edge_1E; Support Region[{GammaC}] ; Entity EdgesOf[All]; }
      { Name se2e ; NameOfCoef we2e; Function BF_Edge_2E; Support Region[{GammaC}] ; Entity EdgesOf[All]; }
    }
  }
}


Formulation {
  { Name Maxwell_e; Type FemEquation;
    Quantity {
      { Name e; Type Local;  NameOfSpace Hcurl_e; }
      { Name ls; Type Local;  NameOfSpace HcurlLS_e; }
      { Name lc; Type Local;  NameOfSpace HcurlLC_e; }
    }
    Equation {
      Galerkin { [ Dof{d e} , {d e} ];
                 In Omega; Integration I1; Jacobian JVol;  }

      Galerkin { [ -k^2 * Dof{e} , {e} ];
                 In Omega; Integration I1; Jacobian JVol;  }

      Galerkin { [ Dof{ls} , {e} ];
                 In GammaS; Integration I1; Jacobian JSur;  }

      Galerkin { [ Dof{lc} , {e} ];
                 In GammaC; Integration I1; Jacobian JSur;  }

      Galerkin { [ Dof{e} , {ls} ];
                 In GammaS; Integration I1; Jacobian JSur;  }

      Galerkin { [ Dof{e} , {lc} ];
                 In GammaC; Integration I1; Jacobian JSur;  }

      Galerkin { [ Vector [0, 1, 0] , {ls} ];
                 In GammaS; Integration I1; Jacobian JSur;  }

      Galerkin { [ Vector [0, 0, 0] , {lc} ];
                 In GammaC; Integration I1; Jacobian JSur;  }
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
      Print[ e, OnElementsOf Omega, File "maxwell.pos"] ;
    }
  }
}
