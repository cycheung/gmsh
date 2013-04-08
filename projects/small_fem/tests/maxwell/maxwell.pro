Group{
  GammaS = Region[5]; // Source
  GammaC = Region[6]; // Conductor -- Reflection
  Omega  = Region[7]; // Omega
}


Function{
  k = 5;
  fs[GammaS] = Vector [0., -1., 0.];
  fc[GammaC] = Vector [0.,  0., 0.];
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

Constraint {
  { Name Dirichlet_e;
   Case {
     { Region GammaC ;
       Type AssignFromResolution;
       NameOfResolution Maxwell_e_DirichletC ; }

     { Region GammaS ;
       Type AssignFromResolution;
       NameOfResolution Maxwell_e_DirichletS ; }
   }
  }
}

FunctionSpace {
  { Name Hcurl_e; Type Form1;
    BasisFunction {
      // Ordre 1 Complet //
      { Name se   ; NameOfCoef ee   ; Function BF_Edge_1E ; Support Region[{Omega,GammaS,GammaC}] ; Entity EdgesOf[All]; }
      { Name se2e ; NameOfCoef we2e ; Function BF_Edge_2E ; Support Region[{Omega,GammaS,GammaC}] ; Entity EdgesOf[All]; }

      // Ordre 2 Complet //
      { Name se3fa; NameOfCoef we3fa; Function BF_Edge_3F_a ; Support Region[{Omega,GammaS,GammaC}] ; Entity FacetsOf[All]; }
      { Name se3fb; NameOfCoef we3fb; Function BF_Edge_3F_b ; Support Region[{Omega,GammaS,GammaC}] ; Entity FacetsOf[All]; }
      { Name se3fc; NameOfCoef we3fc; Function BF_Edge_3F_c ; Support Region[{Omega,GammaS,GammaC}] ; Entity FacetsOf[All]; }

      { Name se4e ; NameOfCoef we4e ; Function BF_Edge_4E   ; Support Region[{Omega,GammaS,GammaC}] ; Entity EdgesOf[All]; }
    }

    Constraint {
      { NameOfCoef ee   ; EntityType EdgesOf  ; NameOfConstraint Dirichlet_e; }
      { NameOfCoef we2e ; EntityType EdgesOf  ; NameOfConstraint Dirichlet_e; }
      { NameOfCoef we3fa; EntityType FacetsOf ; NameOfConstraint Dirichlet_e; }
      { NameOfCoef we3fb; EntityType FacetsOf ; NameOfConstraint Dirichlet_e; }
      { NameOfCoef we3fc; EntityType FacetsOf ; NameOfConstraint Dirichlet_e; }
      { NameOfCoef we4e ; EntityType EdgesOf  ; NameOfConstraint Dirichlet_e; }
    }
  }
}

Formulation {
  { Name Maxwell_e_DirichletC;
    Quantity {
      { Name e; Type Local; NameOfSpace Hcurl_e; }
    }
    Equation {
      Galerkin { [ Dof{e} , {e} ];
	In GammaC; Integration I1; Jacobian JSur;  }
      Galerkin { [ fc[] , {e} ];
	In GammaC; Integration I1; Jacobian JSur;  }
    }
  }

  { Name Maxwell_e_DirichletS;
    Quantity {
      { Name e; Type Local; NameOfSpace Hcurl_e; }
    }
    Equation {
      Galerkin { [ Dof{e} , {e} ];
	In GammaS; Integration I1; Jacobian JSur;  }
      Galerkin { [ fs[] , {e} ];
	In GammaS; Integration I1; Jacobian JSur;  }
    }
  }

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
      { Name A ; NameOfFormulation Maxwell_e ; }
    }
    Operation {
      Generate[A] ; Solve[A] ;
    }
  }

  { Name Maxwell_e_DirichletC;
    System {
      { Name B; NameOfFormulation Maxwell_e_DirichletC; DestinationSystem A; }
    }
    Operation {
      Generate B; Solve B; TransferSolution B;
    }
  }

  { Name Maxwell_e_DirichletS;
    System {
      { Name C; NameOfFormulation Maxwell_e_DirichletS; DestinationSystem A; }
    }
    Operation {
      Generate C; Solve C; TransferSolution C;
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
