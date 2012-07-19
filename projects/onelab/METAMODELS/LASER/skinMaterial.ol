OL.if( OL.get(Parameters/Model/TENEUR) )
  Density = Variable Teneur
    Real MATC "1000/(0.0616*tx(0)+0.938)" ! kg/m3
  Heat Conductivity = Variable Teneur, DensityBis
    Real MATC "tx(1)/1000*(0.454*tx(0)+0.174)" ! W/(mK) 
  Heat Capacity = Variable Teneur
    Real MATC "2500*tx(0)+1700"  ! J/(kg/K)
OL.else
  Density = Real MATC "1000/(0.0616*teneurw+0.938)"   !1022.45 !1048.88 !1200
  Heat Conductivity = Variable DensityBis
    Real MATC "tx(0)/1000*(0.454*teneurw+0.174)" !0.48 !0.30 !0.58
  Heat Capacity = Real MATC "2500*teneurw+1700" !3325.0 !2325.0 !3800.0 
OL.endif