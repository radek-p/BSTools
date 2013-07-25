
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

        /* windows & menus */
#define MENU00           100
#define MENU02           101
#define POPUP00          102
#define POPUP01          103
#define POPUP02          104
#define POPUP03          105
#define KNOTWIN          106
#define CURVEWIN         107

        /* menu 00 widget numbers */
#define intw00DEG        108
#define sw00SELECT       109
#define sw00MOVE         110
#define sw00SCALE        111
#define sw00ROTATE       112
#define sw00SHEAR        113
#define sw00NURBS        114
#define sw00CLOSED       115
#define sw00FUNCTION     116
#define sw00UNIFORM      117
#define btn00REFINE      118
#define btn00RESET       119
#define sw00PANZOOM      120
#define sw00COORD        121
#define sw00MOVEMANY     122

        /* menu 01 widget numbers */
#define sw01CURVE        123
#define sw01POLYLINE     124
#define sw01BEZPOLY      125
#define sw01TICKS        126
#define sw01CONVH        127
#define sw01BASIS        128
#define sw01POLARF       129
#define sw01CURVGR       130
#define sl01CURVGRSC     131
#define intw01CURVGRDENS 132

        /* menu02 widget numbers */
#define btn02FILE        133
#define btn02EDIT        134
#define btn02VIEW        135
#define btn02ABOUT       136

        /* popup 00 widget numbers */
#define btn00pOPEN       137
#define btn00pSAVE       138
#define btn00pSAVEAS     139
#define btn00pEXPORT     140
#define btn00pEXIT       141

        /* popup 01 widget numbers */
#define btn01pOPEN       142
#define btn01pCANCEL     143
#define lb01pDIRLIST     144
#define lb01pFILELIST    145

        /* popup 02 widget numbers */
#define btn02pSAVE       146
#define btn02pCANCEL     147
#define lb02pDIRLIST     148
#define txt02pSAVE_AS    149
#define txted02pFILENAME 150

        /* popup 03 widget numbers */
#define btn03EXIT        151
#define btn03SAVE        152
#define btn03CANCEL      153


extern char txtFile[];
extern char txtOpen[];
extern char txtSave[];
extern char txtSaveAs[];
extern char txtExport[];
extern char txtExit[];
extern char txtCancel[];

extern char b1text[];
extern char b2text[];
extern char b3text[];
extern char b4text[];
extern char sw19text[];
extern char b5text[];
extern char sw0text[];
extern char sw1text[];
extern char sw2text[];
extern char sw3text[];
extern char txtShear[];
extern char intw2text[];
extern char sw4text[];
extern char sw5text[];
extern char sw6text[];
extern char sw7text[];
extern char sw8text[];
extern char sw9text[];
extern char sw10text[];
extern char sw11text[];
extern char sw12text[];
extern char sw13text[];
extern char sw14text[];
extern char intw15text[];
extern char sw16text[];
extern char sw17text[];
extern char sw18text[];

extern char *InfoMsg[];

extern char MsgReconsider[];

extern char ErrMsgRaiseDeg[];
extern char ErrMsgReduceDeg[];
extern char ErrMsgCannotInsert[];
extern char ErrMsgToManyKnots[];
extern char ErrMsgCannotRemove[];
extern char ErrMsgCannotClose[];
extern char ErrMsgCannotSave[];
extern char ErrMsgCannotRead[];
extern char ErrMsgCannotRefine[];

