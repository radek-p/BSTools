
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

        /* windows & menus */
#define CURVEWIN         100
#define MENU00           101
#define MENU01           102
#define KNOTWIN          103
#define POPUP00          104
#define POPUP01          105
#define POPUP02          106
#define POPUP03          107

        /* menu 00 widget numbers */
#define btn00FILE        108
#define btn00EDIT        109
#define btn00VIEW        110
#define btn00ABOUT       111

        /* menu 01a widget numbers */
#define intw01aDEGREE    112
#define sw01aSELECT      113
#define sw01aMOVE        114
#define sw01aSCALE       115
#define sw01aROTATE      116
#define sw01aSHEAR       117
#define sw01aNURBS       118
#define sw01aCLOSED      119
#define sw01aUNIFORM     120
#define btn01aREFINE     121
#define btn01aRESET      122
#define sw01aPANZOOM     123
#define sw01aCOORD       124
#define sw01aMOVEMANY    125

        /* menu 01b widget numbers */
#define sw01bCPOLY       126
#define sw01bCURVE       127
#define sw01bTICKS       128
#define sw01bBEZPOLY     129
#define sw01bCONVH       130
#define sw01bDIAGF       131
#define sw01bCURVGRAPH   132
#define sl01bCURVSCALE   133
#define sw01bTORSGRAPH   134
#define sl01bTORSSCALE   135
#define intw01bGRAPHDENS 136

        /* popup 00 widget numbers */
#define btn00pOPEN       137
#define btn00pSAVE       138
#define btn00pSAVEAS     139
#define btn00pEXIT       140

        /* popup 01 widget numbers */
#define btn01pOPEN       141
#define btn01pCANCEL     142
#define lb01pDIRLIST     143
#define lb01pFILELIST    144

        /* popup 02 widget numbers */
#define btn02pSAVE       145
#define btn02pCANCEL     146
#define lb02pDIRLIST     147
#define txt02pSAVE_AS    148
#define txted02pFILENAME 149

        /* popup 03 widget numbers */
#define btn03EXIT        150
#define btn03SAVE        151
#define btn03CANCEL      152


extern char txtFile[];
extern char txtOpen[];
extern char txtSave[];
extern char txtSaveAs[];
extern char txtExit[];
extern char txtCancel[];

extern char txtEdit[];
extern char txtView[];
extern char txtAbout[];
extern char txtUniform[];
extern char txtRefine[];
extern char txtReset[];

extern char txtDegree[];
extern char txtSelectUnselect[];
extern char txtMove[];
extern char txtScale[];
extern char txtRotate[];
extern char txtShear[];
extern char txtControlPolygon[];
extern char txtCurve[];
extern char txtTicks[];
extern char txtBezierPolygons[];
extern char txtConvexHulls[];
extern char txtNURBS[];
extern char txtClosed[];
extern char txtDiagonalForms[];
extern char txtCurvatureGraph[];
extern char txtTorsionGraph[];
extern char txtGraphDensity[];
extern char txtPanZoom[];
extern char txtCoordinates[];
extern char txtMoveManyKnots[];

extern char *InfoMsg[5];

extern char MsgReconsider[];

extern char ErrMsgCannotSave[];
extern char ErrMsgCannotRead[];
extern char ErrMsgCannotRefine[];

