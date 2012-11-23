// Switched Reluctance Motor Parameter File (2D)

DefineConstant[Flag_Symmetry = {1,Choices {0,1}}] ;//0 => Full model; 1 => Half model

//------------------------------------------------------------
//------------------------------------------------------------

Ns = 6.;  // Num. of Stator Poles
Nr = 4.;  // Num. of Rotor Poles
AxialLength = 60e-3 ;

//------------------------------------------------------------
// Physical regions
//------------------------------------------------------------

STATOR = 1000 ;
ROTOR  = 2000 ;
SHAFT  = 2100 ;

COILN  = 1100 ;
COILP  = 1200 ;

BND_STATOR = 5555 ;

AIRGAP_STATOR  = 3000 ;
AIRGAP_ROTOR   = 3100 ;
AIRGAP_ROTORIN = 3101 ;
AIRGAP_DUMMY   = 3102 ;

LMB_STATOR  = 1111 ;
LMB_ROTOR   = 2222 ;
LMB_ROTORAUX = 2223 ;
PMB_ROTORREF = 222 ;
PMB_ROTORAUX = 223 ;

LSTATOR_RIGHT = 3333 ;
LSTATOR_LEFT  = 3334 ;
LROTOR_RIGHT  = 4444 ;
LROTOR_LEFT   = 4445 ;

LINE_NICEVIEW = 77777 ;
