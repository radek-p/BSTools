
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define win3D0                    100 /* a multiple of 4 */
#define MENU0                     102
#define MENU1                     103
#define MENU2                     104
#define POPUP00                   105
#define POPUP01                   106
#define POPUP02                   107
#define POPUP03                   108
#define POPUP04                   109
#define STATUSLINE0               110
#define sw02STATUS                111

#define win2D0                    112
#define MENU3                     113
#define MENU4                     114
#define MENU5                     115
#define MENU6                     116
#define MENU7                     117
#define MENU8                     118
#define MENU9                     119
#define MENU10                    120
#define MENU11                    121
#define MENU12                    122
#define MENU13                    123
#define MENU14                    124
#define STATUSLINE1               125
#define sw16STATUS                126
#define POPUP10                   127
#define POPUP11                   128
#define POPUP12                   129

/* widget numbers, window 0, menu 0 */
#define btn00aRESET               130
#define sw00aMARK                 131
#define sw00aMOVE                 132
#define sw00aSCALE                133
#define sw00aROTATE               134
#define sw00aSHEAR                135
#define sw00aPAN_ZOOM             136
#define sw00aCOORDINATES          137

#define sw00bCONTROLNET           138
#define sw00bSURFACE              139
#define sw00bBEZIER_NETS          140
#define intw00bU_DENSITY          141
#define intw00bV_DENSITY          142
#define sw00bCONSTRPOLY           143
#define sw00bTRANSFNET            144
#define sw00bPAN_ZOOM             145
#define sw00bCOORDINATES          146

#define btn00cRENDER_STOP         147
#define sw00cGAUSSIAN_C           148
#define sw00cMEAN_C               149
#define sw00cLAMBERTISO_C         150
#define sw00cREFLECTION_C         151
#define sw00cHIGHLIGHT_C          152
#define sw00cSECTIONS_C           153
#define sw00cPARAM_C              154
#define sw00cGAUSSIAN_D           155
#define sw00cMEAN_D               156
#define sw00cLAMBERTISO_D         157
#define sw00cREFLECTION_D         158
#define sw00cHIGHLIGHT_D          159
#define sw00cSECTIONS_D           160
#define sw00cPARAM_D              161
#define sl00cREND_CFRANGE         162
#define sl00cREND_DFSF            163
#define sw00cSHADOWS              164
#define sw00cANTIALIAS            165

#define btn00dRENDER_STOP         166
#define sw00dLIGHT0DIR            167
#define sl00dLIGHT0INT            168
#define sw00dLIGHT1DIR            169
#define sl00dLIGHT1INT            170
#define sw00dLIGHT2DIR            171
#define sl00dLIGHT2INT            172
#define sw00dLIGHT3DIR            173
#define sl00dLIGHT3INT            174
#define sl00dLIGHTAMB             175
#define sw00dREFLECTIONFRAME      176
#define sw00dHIGHLIGHTFRAME       177
#define sw00dSECTIONSFRAME        178

/* widget numbers, window 0, menu 1 */
#define btn01FILEpopup            179
#define btn01EDITpopup            180
#define btn01VIEW                 181
#define btn01PICTURE              182
#define btn01ABOUT                183

/* widget numbers, window 0, menu 2 (popup) */
#define btn02OPENpopup            184
#define btn02SAVEpopup            185
#define btn02SAVEASpopup          186
#define btn02EXPORT               187
#define btn02EXIT                 188

/* widget numbers, window 0, menu 3 (popup, Open) */
#define txt03DIRECTORY            189
#define btn03OPEN                 190
#define btn03CANCEL               191
#define lb03DIRLIST               192
#define lb03FILELIST              193

/* widget numbers, window 0, menu 4 (popup, Save as) */
#define txt04DIRECTORY            194
#define btn04SAVE                 195
#define btn04CANCEL               196
#define lb04DIRLIST               197
#define txt04SAVE_AS              198
#define txted04FILENAME           199

/* widget numbers, window 0 menu 5 (popup, Edit) */
#define btn05SURFACE              200
#define btn05LIGHT                201

/* widget numbers, window 0 menu 6 (popup, Exit) */
#define btn06EXIT                 202
#define btn06SAVE                 203
#define btn06CANCEL               204

/* widget numbers, window 1, menu 3 */
#define txt13GENERAL              205
#define intw13DEGREE_U            206
#define intw13DEGREE_V            207
#define sw13CLOSED_U              208
#define sw13CLOSED_V              209
#define btn13EQUIDIST_U           210
#define btn13EQUIDIST_V           211
#define btn13FLIP                 212
#define sw13DOMAIN_NET            213
#define sw13MARK_UNMARK           214
#define sw13MOVE_MANY_KNOTS       215
#define sw13PAN_ZOOM              216
#define sw13COORDINATES           217

/* widget numbers, window 1, menu 4 */
#define btn14aTYPE                218

#define btn14bTYPE                219
#define btn14bEDIT                220
#define btn14bVIEW                221
#define btn14bARCS                222

/* widget numbers, window 1, menu 5 */
#define txt15aSPHERICAL           223
#define sw15aBIND                 224
#define sw15aEQUATOR              225
#define sw15aMERIDIAN             226
#define btn15aRESET               227
#define intw15aDEGREE             228
#define sw15aCLOSED               229
#define sw15aNURBS                230
#define sw15aMARK_UNMARK          231
#define sw15aMOVE                 232
#define sw15aSCALE                233
#define sw15aROTATE               234
#define sw15aSHEAR                235
#define sw15aPAN_ZOOM             236
#define sw15aCOORDINATES          237
#define sw15aMOVE_MANY_KNOTS      238

#define sw15bCONTROL_POLYGON      239
#define sw15bCURVE                240
#define sw15bBEZIER_POLYGONS      241
#define sw15bTICKS                242

#define btn15cQUARTER_CIRCLE      243
#define btn15cHALF_CIRCLE         244
#define btn15cFULL_CIRCLE         245
#define dial15cARC_ANGLE          246

/* widget numbers, window 1, menu 6 (popup) */
#define btn16GENERAL              247
#define btn16SPHERICAL            248
#define btn16SWEPT                249
#define btn16BLENDING             250

/* widget numbers, window 1, menu 9 - swept surface */
/* widget numbers, window 1, menu 10 - swept surface */
#define txt110SWEPT               251

/* widget numbers, window 1, menu 11 - blending surface */
#define btn111TRANSFORM           252
#define btn111OPTIONS             253
#define btn111INFO                254

/* widget numbers, window 1, menu 12 - blending surface */
#define txt112BLENDING            255
#define sw112BLENDING_G1          256
#define sw112BLENDING_G2          257
#define btn112INIT_BLENDING       258
#define btn112BLENDING_REFINE     259
#define sw112CLAMPED_BLENDING     260
#define sw112CLOSED_BLENDING      261
#define sw112CONSTRAINTS          262
#define sw112TRIHARMONIC_BLENDING 263
#define btn112BLENDING_LMT_ITER   264
#define sl112NONLIN_BLENDING_C    265
#define sw112BLENDING_OPT_ENTIRE  266
#define intw112BLENDING_UMIN      267
#define intw112BLENDING_UMAX      268
#define intw112BLENDING_VMIN      269
#define intw112BLENDING_VMAX      270

/* widget numbers, window1, popup11 - pretransformation for blending surfaces */
#define txtedP11A0                271
#define txtedP11A1                272
#define txtedP11A2                273
#define txtedP11A3                274
#define txtedP11A4                275
#define txtedP11A5                276
#define txtedP11A6                277
#define txtedP11A7                278
#define txtedP11A8                279
#define btnP11IDENTITY            280
#define btnP11OK                  281

/* widget numbers, window1, popup12 - blending optimization options */
#define intwP12BLENDING_LMT_NITER 282
#define intwP12BLENDING_QUAD1     283
#define intwP12BLENDING_QUAD2     284
#define swP12SHOWSTEPS            285
#define swP12DUMPDATA             286
#define btnP12OK                  287


/* widget description character strings; these might be translated */
extern char txtNULL[];
extern char txtReset[];
extern char txtInit[];
extern char txtEdit[];
extern char txtView[];
extern char txtPicture[];
extern char txtAbout[];
extern char txtDegree[];
extern char txtDegree_u[];
extern char txtDegree_v[];
extern char txtEquidist_u[];
extern char txtEquidist_v[];
extern char txtFlip[];
extern char txtPan_zoom[];
extern char txtCoordinates[];
extern char txtControl_net[];
extern char txtSurface[];
extern char txt_Surface[];
extern char txt_Light[];
extern char txtMove[];
extern char txtScale[];
extern char txtRotate[];
extern char txtShear[];
extern char txtU_density[];
extern char txtV_density[];
extern char txtMark_unmark[];
extern char txtBezier_nets[];
extern char txtConstrPoly[];
extern char txtTransfNet[];
extern char txtMove_many_knots[];
extern char txtClosed[];
extern char txtClosed_u[];
extern char txtClosed_v[];
extern char txtDomain_net[];
extern char txtFile[];
extern char txtOpen[];
extern char txtSave[];
extern char txtSave_as[];
extern char txtExport[];
extern char txtExit[];
extern char txtOK[];
extern char txtCancel[];
extern char txtSurface_Type[];
extern char txtGeneral[];
extern char txtGeneral_BSpline[];
extern char txtSpherical[];
extern char txtSpherical_product[];
extern char txtSwept[];
extern char txtBlending[];
extern char txtEquator[];
extern char txtMeridian[];
extern char txtBind[];
extern char txtNURBS[];
extern char txtCurve[];
extern char txtControl_polygon[];
extern char txtBezier_polygons[];
extern char txtTicks[];
extern char txtRender[];
extern char txtC[];
extern char txtD_shape_func[];
extern char txtGaussian[];
extern char txtMean[];
extern char txtIsophotes[];
extern char txtReflection[];
extern char txtHighlight[];
extern char txtSections[];
extern char txtParamQual[];
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
extern char txtArcs[];
extern char txtQuarter[];
extern char txtHalf[];
extern char txtFull_circle[];
extern char txtArc_angle[];
extern char txtSwept_Surface[];
extern char txtBlending_Surface[];
extern char txtClamped[];
extern char txtBiharmonic[];
extern char txtTriharmonic[];
extern char txtRefine[];
extern char txtIterate[];
extern char txtConstraints[];
extern char txtTransform[];
extern char txtOptions[];
extern char txtInfo[];
extern char txtOptimize[];
extern char txtInterrupt[];
extern char txtIterationLimit[];
extern char txtQuadratureKnots[];
extern char txtShowSteps[];
extern char txtDumpData[];
extern char txtIdentity[];
extern char txtG1[];
extern char txtG2[];

/* message texts */
extern char *InfoMsg[];

extern char MsgReconsider[];

extern char ErrorMsgCannotOpen[];
extern char ErrorMsgCannotSave[];
extern char ErrorMsgCannotRaiseDeg[];
extern char ErrorMsgCannotReduceDeg[];
extern char ErrorMsgCannotInsertKnot[];
extern char ErrorMsgTooManyKnots[];
extern char ErrorMsgCannotRemoveKnot[];
extern char ErrorMsgNotEnoughKnots[];
extern char ErrorMsgNotEnoughUKnots[];
extern char ErrorMsgNotEnoughVKnots[];
extern char ErrorMsgFileWritingError[];
extern char ErrorMsgRefineBlending[];
extern char ErrorMsgBindBlending[];
extern char ErrorMsgNonlinBlending[];
extern char ErrorMsgCannotAddConstraint[];
extern char ErrorMsgIncorrectNumber[];
extern char ErrorMsgBiquadratic[];
extern char ErrorMsgBicubic[];

extern char ErrorChildProcessNotActive[];
extern char ErrorChildProcessTerminated[];


extern char ErrorMsgNotImplemented[];

extern char WarningMsgPovExport[];

