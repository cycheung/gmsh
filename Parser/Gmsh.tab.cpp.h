typedef union {
  char    *c;
  int      i;
  unsigned int u;
  double   d;
  double   v[5];
  Shape    s;
  List_T  *l;
} YYSTYPE;
#define	tDOUBLE	258
#define	tSTRING	259
#define	tBIGSTR	260
#define	tEND	261
#define	tAFFECT	262
#define	tDOTS	263
#define	tPi	264
#define	tExp	265
#define	tLog	266
#define	tLog10	267
#define	tSqrt	268
#define	tSin	269
#define	tAsin	270
#define	tCos	271
#define	tAcos	272
#define	tTan	273
#define	tAtan	274
#define	tAtan2	275
#define	tSinh	276
#define	tCosh	277
#define	tTanh	278
#define	tFabs	279
#define	tFloor	280
#define	tCeil	281
#define	tFmod	282
#define	tModulo	283
#define	tHypot	284
#define	tPrintf	285
#define	tPoint	286
#define	tCircle	287
#define	tEllipsis	288
#define	tLine	289
#define	tSurface	290
#define	tSpline	291
#define	tVolume	292
#define	tCharacteristic	293
#define	tLength	294
#define	tParametric	295
#define	tElliptic	296
#define	tPlane	297
#define	tRuled	298
#define	tTransfinite	299
#define	tComplex	300
#define	tPhysical	301
#define	tUsing	302
#define	tBump	303
#define	tProgression	304
#define	tRotate	305
#define	tTranslate	306
#define	tSymmetry	307
#define	tDilate	308
#define	tExtrude	309
#define	tDuplicata	310
#define	tLoop	311
#define	tInclude	312
#define	tRecombine	313
#define	tDelete	314
#define	tCoherence	315
#define	tView	316
#define	tAttractor	317
#define	tLayers	318
#define	tScalarTetrahedron	319
#define	tVectorTetrahedron	320
#define	tTensorTetrahedron	321
#define	tScalarTriangle	322
#define	tVectorTriangle	323
#define	tTensorTriangle	324
#define	tScalarLine	325
#define	tVectorLine	326
#define	tTensorLine	327
#define	tScalarPoint	328
#define	tVectorPoint	329
#define	tTensorPoint	330
#define	tBSpline	331
#define	tNurbs	332
#define	tOrder	333
#define	tWith	334
#define	tBounds	335
#define	tKnots	336
#define	tColor	337
#define	tFor	338
#define	tEndFor	339
#define	tScript	340
#define	tExit	341
#define	tMerge	342
#define	tB_SPLINE_SURFACE_WITH_KNOTS	343
#define	tB_SPLINE_CURVE_WITH_KNOTS	344
#define	tCARTESIAN_POINT	345
#define	tTRUE	346
#define	tFALSE	347
#define	tUNSPECIFIED	348
#define	tU	349
#define	tV	350
#define	tEDGE_CURVE	351
#define	tVERTEX_POINT	352
#define	tORIENTED_EDGE	353
#define	tPLANE	354
#define	tFACE_OUTER_BOUND	355
#define	tEDGE_LOOP	356
#define	tADVANCED_FACE	357
#define	tVECTOR	358
#define	tDIRECTION	359
#define	tAXIS2_PLACEMENT_3D	360
#define	tISO	361
#define	tENDISO	362
#define	tENDSEC	363
#define	tDATA	364
#define	tHEADER	365
#define	tFILE_DESCRIPTION	366
#define	tFILE_SCHEMA	367
#define	tFILE_NAME	368
#define	tMANIFOLD_SOLID_BREP	369
#define	tCLOSED_SHELL	370
#define	tADVANCED_BREP_SHAPE_REPRESENTATION	371
#define	tFACE_BOUND	372
#define	tCYLINDRICAL_SURFACE	373
#define	tCONICAL_SURFACE	374
#define	tCIRCLE	375
#define	tTRIMMED_CURVE	376
#define	tGEOMETRIC_SET	377
#define	tCOMPOSITE_CURVE_SEGMENT	378
#define	tCONTINUOUS	379
#define	tCOMPOSITE_CURVE	380
#define	tTOROIDAL_SURFACE	381
#define	tPRODUCT_DEFINITION	382
#define	tPRODUCT_DEFINITION_SHAPE	383
#define	tSHAPE_DEFINITION_REPRESENTATION	384
#define	tELLIPSE	385
#define	tTrimmed	386
#define	tSolid	387
#define	tEndSolid	388
#define	tVertex	389
#define	tFacet	390
#define	tNormal	391
#define	tOuter	392
#define	tLoopSTL	393
#define	tEndLoop	394
#define	tEndFacet	395
#define	tAND	396
#define	tOR	397
#define	tNOTEQUAL	398
#define	tEQUAL	399
#define	tAPPROXEQUAL	400
#define	tAFFECTPLUS	401
#define	tAFFECTMINUS	402
#define	tAFFECTTIMES	403
#define	tAFFECTDIVIDE	404
#define	tLESSOREQUAL	405
#define	tGREATEROREQUAL	406
#define	tCROSSPRODUCT	407
#define	UNARYPREC	408
#define	tPLUSPLUS	409
#define	tMINUSMINUS	410


extern YYSTYPE yylval;
