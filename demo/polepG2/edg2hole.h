
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#define MAIN0_RECT_NUM   4
#define MENU0A_RECT_NUM 21
#define MENU0B_RECT_NUM  7

#define MAIN1_RECT_NUM   7
#define MENU1A_RECT_NUM 12
#define MENU1B_RECT_NUM 39
#define MENU1C_RECT_NUM 15

#define STATE_MOVINGPOINT    3
#define STATE_MOVINGKNOT     4
#define STATE_TURNINGVIEWER  5
#define STATE_ZOOM           6
#define STATE_PARZOOM        7
#define STATE_RESIZE_X       8
#define STATE_RESIZE_Y       9
#define STATE_RESIZE_XY     10
#define STATE_MOVING_DOMCP  11
#define STATE_DOMZOOM       12
#define STATE_MOVING_FIRSTP 13
#define STATE_THINKING0     14
#define STATE_THINKING1     15
#define STATE_MOVING_CENTP  16

#define FIRST_GEOMW 1

#define MAX_ZOOM    100.0
#define MIN_ZOOM      1.0
#define MAX_PARZOOM  10.0
#define MIN_PARZOOM   0.1
#define MAX_DOMZOOM 100.0
#define MIN_DOMZOOM   1.0

extern int main_rect_num;
extern int menu0, menu1;

extern ed_rect edrect0[MAIN0_RECT_NUM];
extern ed_rect menurect0a[MENU0A_RECT_NUM];
extern ed_rect menurect0b[MENU0B_RECT_NUM];

extern ed_rect edrect1[MAIN1_RECT_NUM];
extern ed_rect menurect1a[MENU1A_RECT_NUM];
extern ed_rect menurect1b[MENU1B_RECT_NUM];
extern ed_rect menurect1c[MENU1C_RECT_NUM];

extern int hole_k;
extern boolean HoleKSwitch[4];

extern boolean swDisplayFirstPartition, swDisplayDomainCP,
  swDisplayDomSurrPatches, swDisplayDomCurves, swDisplayDomAuxPatches,
  swDisplayDomPatches, swDisplayPartition, swDisplayDomNumbers;

extern boolean swDisplayCentralPoint, swUseDerivatives1, swAltDomCurves,
  swRestrictBasis;

extern boolean swDisplaySurfCP, swDisplaySurfPatches, swDisplayFinalPatches,
  swDisplaySurfNumbers, swDisplayNLFinalPatches;
extern boolean PictureIsOn;

void SetHoleK ( int k );
void SetConstraintSwitches ( int k );
void SetDomainParam ( void );
void DomWinMsg ( ed_rect *er, int msg, int key, int x, int y );
void DrawDomWin ( ed_rect *er );
void KnotWinMsg ( ed_rect *er, int msg, int key, int x, int y );
void DrawKnotWin ( ed_rect *er );
void ResizeWindow0 ( void );
void FixDomainCP ( void );

void SetSurfParam ( void );
void TurnFinalSurface ( void );
void TurnNLFinalSurface ( void );
void CompGeomWinSizes ( void );
void RedrawGeomWindows ( void );
void DrawParWindow ( ed_rect *er );
void ParWindowMsg ( ed_rect *er, int msg, int key, int x, int y );
void DrawPerspWindow ( ed_rect *er );
void SetZoom ( float zf );
void PerspWindowMsg ( ed_rect *er, int msg, int key, int x, int y );
void DrawSpecialWindow ( ed_rect *er );
void SpecialMsg ( ed_rect *er, int msg, int key, int x, int y );
void ThinkingMsg ( ed_rect *er, int msg, int key, int x, int y );
void ResizeWindow1 ( boolean reset );
void FitToDomain ( void );
void Flatten ( void );
void FitToSurface ( void );
void OnRendering ( void );
void OffRendering ( void );
void OnOffRendering ( void );
void SwitchMenu0 ( void );
void SwitchMenu1 ( int menu );

void WriteFile ( void );

