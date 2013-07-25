
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#define win3D0                100
#define MENU0                 102
#define MENU1                 103
#define MENU2                 104
#define POPUP00               105
#define POPUP01               106
#define POPUP02               107
#define POPUP03               108
#define POPUP04               109
#define STATUSLINE0           110
#define sw02STATUS            111

#define win2D0                112
#define win1D0                113
#define MENU3                 114
#define MENU4                 115
#define MENU5                 116
#define STATUSLINE1           117
#define sw15STATUS            118

#define btn00FILE             119
#define btn00DATA             120
#define btn00EDIT             121
#define btn00VIEW             122
#define btn00PICTURE          123
#define btn00ABOUT            124

#define intw01aHOLE_SIDES     125
#define sl01aPARAM0           126
#define sl01aPARAM1           127
#define sl01aPARAM2           128
#define sl01aPARAM3           129
#define sl01aPARAM4           130

#define sw01bMARK_UNMARK      131
#define sw01bMOVE             132
#define sw01bSCALE            133
#define sw01bROTATE           134
#define sw01bSHEAR            135
#define sw01PAN_ZOOM          136
#define sw01COORDINATES       137

#define sw01cMARK_UNMARK      138
#define sw01cMOVE             139
#define sw01cSCALE            140
#define sw01cROTATE           141
#define sw01cSHEAR            142
#define sw01c1ST_SURF_CONSTR  143
#define sw01c2ND_SURF_CONSTR  144
#define sw01c1ST_CONSTR       145
#define sw01c2ND_CONSTR       146
#define sw01cZERO_CONSTR      147
#define btn01cGET_CURRENT     148
#define scw01cCONSTR_SWITCHES 149
#define menu01cSCROLLED_SW    150
#define sw01cCONSTR_SWITCH    151
#define MAX_CONSTR_SWITCH_ROWS     16
#define MAX_CONSTR_SWITCH_COLS     16
#define NUM_CONSTR_SWITCHES       257  /* 16*16+1 */

#define sw01dCPOINTS          408 /* 150+257 */
#define sw01dSURFACE          409
#define sw01dFIRST            410
#define sw01dSECOND           411
#define sw01dNUMBERS          412
#define sw01dCONSTRAINT_FRAME 413

#define btn01eRENDER_STOP     414
#define sw01eGAUSSIAN_C       415
#define sw01eMEAN_C           416
#define sw01eLAMBERTISO_C     417
#define sw01eREFLECTION_C     418
#define sw01eHIGHLIGHT_C      419
#define sw01eSECTIONS_C       420
#define sw01eGAUSSIAN_D       421
#define sw01eMEAN_D           422
#define sw01eLAMBERTISO_D     423
#define sw01eREFLECTION_D     424
#define sw01eHIGHLIGHT_D      425
#define sw01eSECTIONS_D       426
#define sw01eSHADOWS          427
#define sw01eANTIALIAS        428

#define btn01fRENDER_STOP     429
#define sw01fLIGHT0DIR        430
#define sl01fLIGHT0INT        431
#define sw01fLIGHT1DIR        432
#define sl01fLIGHT1INT        433
#define sw01fLIGHT2DIR        434
#define sl01fLIGHT2INT        435
#define sw01fLIGHT3DIR        436
#define sl01fLIGHT3INT        437
#define sl01fLIGHTAMB         438
#define sw01fREFLECTIONFRAME  439
#define sw01fHIGHLIGHTFRAME   440
#define sw01fSECTIONSFRAME    441

#define btnP00NEW             442
#define btnP00OPEN            443
#define btnP00SAVE            444
#define btnP00SAVEAS          445
#define btnP00EXIT            446

#define txtP01DIRNAME         447
#define btnP01OPEN            448
#define btnP01CANCEL          449
#define lbP01DIRLIST          450
#define lbP01FILELIST         451

#define txtP02DIRNAME         452
#define btnP02SAVE            453
#define btnP02CANCEL          454
#define lbP02DIRLIST          455
#define txtP02SAVE_AS         456
#define txtedP02FILENAME      457

#define btnP03SURFACE         458
#define btnP03CONSTRAINTS     459
#define btnP03LIGHT           460

#define btnP04EXIT            461
#define btnP04SAVE            462
#define btnP04CANCEL          463

#define btn10OPTIONS          464
#define btn10DATA             465
#define btn10EDIT             466
#define btn10VIEW             467
#define btn10INFO             468

#define sw11a1FIRST           469
#define sw11a1SECOND          470
#define intw11a1ORDER         471
#define sw11a1RESTRICTED      472
#define sw11a1COONS           473
#define sw11a1BEZIER          474
#define sw11a1SPLINE          475
#define intw11a1NK            476
#define intw11a1M1            477
#define intw11a1M2            478
#define sw11a1LINEAR          479
#define sw11a1QUASIG2         480
#define sl11a1QUASIG2CONST    481
#define sw11a1ALTCENTRE       482
#define sw11a1GAUSSLEGENDRE   483

#define sw11a2FIRST           484
#define sw11a2SECOND          485
#define intw11a2ORDER         486
#define sw11a2RESTRICTED      487
#define sw11a2COONS           488
#define sw11a2BEZIER          489
#define sw11a2SPLINE          490
#define intw11a2NK            491
#define intw11a2M1            492
#define intw11a2M2            493
#define sw11a2LINEAR          494
#define sw11a2QUASIG2         495
#define sl11a2QUASIG2CONST    496
#define sw11a2ALTCENTRE       497
#define sw11a2GAUSSLEGENDRE   498

#define intw11bHOLE_SIDES     499
#define sl11bPARAM0           500
#define sl11bPARAM1           501
#define sl11bPARAM2           502
#define sl11bPARAM3           503

#define sw11cMARK_UNMARK      504
#define sw11cMOVE             505
#define sw11cSCALE            506
#define sw11cROTATE           507
#define sw11cSHEAR            508
#define sw11PAN_ZOOM          509
#define sw11COORDINATES       510

#define sw11dDOMCPOINTS       511
#define sw11dDOMSURRPATCHES   512
#define sw11dDOMPATCHES1      513
#define sw11dDOMPATCHES2      514
#define sw11dDOMNUMBERS       515


extern char txtNULL[];
extern char txtFile[];
extern char txtData[];
extern char txtEdit[];
extern char txtConstraints[];
extern char txtLight[];
extern char txtView[];
extern char txtPicture[];
extern char txtAbout[];
extern char txtSides[];
extern char txtMark_unmark[];
extern char txtMove[];
extern char txtScale[];
extern char txtRotate[];
extern char txtShear[];
extern char txtPan_zoom[];
extern char txtCoordinates[];
extern char txtNew[];
extern char txtOpen[];
extern char txtSave[];
extern char txtSave_as[];
extern char txtCancel[];
extern char txtExit[];
extern char txtReset[];
extern char txtOptions[];
extern char txtInfo[];
extern char txtControl_net[];
extern char txtSurr_patches[];
extern char txtDomain_patches1[];
extern char txtDomain_patches2[];
extern char txtNumbers[];
extern char txtFirst[];
extern char txtSecond[];
extern char txtConstrFrame[];
extern char txtGCOrder[];
extern char txtRestricted[];
extern char txtCoons[];
extern char txtBezier[];
extern char txtSpline[];
extern char txtNknots[];
extern char txtMult1[];
extern char txtMult2[];
extern char txtSurface[];
extern char txt_surface[];
extern char txtLinear[];
extern char txtQuasiG2[];
extern char txtAltCentre[];
extern char txtGaussLegendre[];
extern char txtRender[];
extern char txtStop[];
extern char txtC[];
extern char txtD_shape_func[];
extern char txtGaussian[];
extern char txtMean[];
extern char txtIsophotes[];
extern char txtReflection[];
extern char txtHighlight[];
extern char txtSections[];
extern char txtShadows[];
extern char txtAntialias[];
extern char txtLight_0[];
extern char txtLight_1[];
extern char txtLight_2[];
extern char txtLight_3[];
extern char txtAmbient[];
extern char txtReflectionFrame[];
extern char txtHighlightFrame[];
extern char txtSectionsFrame[];
extern char txtFirstType[];
extern char txtSecondType[];
extern char txtZero[];
extern char txtCurrent[];
extern char txtCentralPoint[];

extern char *InfoMsg[];

extern char MsgReconsider[];

extern char ErrMsgNoAction[];
extern char ErrMsgBadFilename[];
extern char ErrMagBadConstrFrame[];
extern char ErrMsgFileReading[];
extern char ErrMsgFileWriting[];

