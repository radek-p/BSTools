
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005,2007                             */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "raybez.h"
#include "eg1holef.h"

#include "oldxgedit.h"
#include "datagenf.h"
#include "g1ekernel.h"
#include "edg1hole.h"
#include "render.h"

int main_rect_num = MAIN0_RECT_NUM;

ed_rect edrect0[MAIN0_RECT_NUM];
ed_rect menurect0a[MENU0A_RECT_NUM];
ed_rect menurect0b[MENU0B_RECT_NUM];

ed_rect edrect1[MAIN1_RECT_NUM];
ed_rect menurect1a[MENU1A_RECT_NUM];
ed_rect menurect1b[MENU1B_RECT_NUM];
ed_rect menurect1c[MENU1C_RECT_NUM];

int menu0 = 0;  /* which menu is shown? */
int menu1 = 0;


void GetTextTitle ( ed_rect *er, char **title )
{
  static char s[2];

  switch ( CurrentWindow () ) {
case 0:
    switch ( er->id ) {
  case 5: sprintf ( s, "%1d", hole_k ); *title = s;  break;
  default: break;
    }
    break;

case 1:
    switch ( er->id ) {
  case 30: *title = "Curvature";  break;
  case 33: *title = "Lines";      break;
  case 40: *title = "Edit";       break;
  default: *title = "";           break;
    }
    break;

default:
    break;
  }
} /*GetTextTitle*/

void DisplayHoleK ( ed_rect *er )
{
  char s[2];

  XSetForeground ( thedisplay, thegc, c_green );
  XFillRectangle ( thedisplay, thepixmap, thegc, er->x, er->y, er->w, er->h );
  sprintf ( s, "%1d", hole_k );
  XSetForeground ( thedisplay, thegc, c_white );
  XDrawString ( thedisplay, thepixmap, thegc, er->x+8, er->y+14, s, strlen(s) );
  XCopyArea ( thedisplay, thepixmap, thewindow, thegc,
              er->x, er->y, er->w, er->h, er->x, er->y );
} /*DisplayHoleK*/

void SetConstraintSwitches ( int k )
{
  int i;

  swConstraintsOn = swZeroDer = swNormalConstr = swPointMode = false;
  for ( i = 0; i < 33; i++ )
    swConstraint[i] = swZConstraint[i] = false;
  constrql1 = constrql2 = 0;
  for ( i = 1; i < 1+4*k; i++ ) {
    if ( (i-1) % 4 < 2 ) {
      menurect1b[i+1].msgproc = SwitchMsg;
      menurect1b[i+1].redraw = DrawSwitch;
    }
    else {
      menurect1b[i+1].msgproc = EmptyMsg;
      menurect1b[i+1].redraw = DrawEmpty;
    }
  }
  for ( i = 1+4*k; i < 33; i++ ) {
    menurect1b[i+1].msgproc = EmptyMsg;
    menurect1b[i+1].redraw = DrawEmpty;
  }
} /*SetConstraintSwitches*/

void MySwitchProc ( ed_rect *er )
{
  switch ( CurrentWindow () ) {
case 0:
    switch ( er->id ) {
  case 1: SetHoleK ( 3 );  goto updslidebars;
  case 2: SetHoleK ( 5 );  goto updslidebars;
  case 3: SetHoleK ( 6 );  goto updslidebars;
  case 4: SetHoleK ( 8 );
updslidebars:
    if ( hole_k == 6 ) {
      menurect1a[3].msgproc = SlidebarMsg;
      menurect1a[3].redraw  = DrawSlidebar;
    }
    else {
      menurect1a[3].msgproc = EmptyMsg;
      menurect1a[3].redraw  = DrawEmpty;
    }
    PictureIsOn = false;
    redraw_all ();
    break;

  case 8:  case 9:  case 10:  case 11:  case 12:  case 13:  case 14:
    edrect0[0].redraw ( &edrect0[0] );
    er->redraw ( er );
    break;

  case 18:
    if ( swDisplayCentralPoint ) {
      redraw ();
      break;
    }
    else {
      swUseDerivatives1 = swAltDomCurves = false;
      goto cont2;
    }

  case 19:
    if ( swUseDerivatives1 ) {
      swDisplayCentralPoint = true;
      redraw ();
      break;
    }
    else goto cont2;

  case 21:
    if ( swAltDomCurves )
      swDisplayCentralPoint = true;
    goto cont2;

  case 22:
cont2:
    swDisplayFinalPatches = false;
    RecreateDomain ();
    PictureIsOn = false;
    redraw_all ();
    break;

  case 23:
    swBezierPatches = !swCoonsPatches;
    goto switchspace;

  case 24:
    swCoonsPatches = !swBezierPatches;
switchspace:
    edrect0[2].redraw ( &edrect0[2] );
    FinalSurfValid = PictureIsOn = false;
    if ( swDisplayNLFinalPatches && !swDisplayFinalPatches ) {
      swDisplayNLFinalPatches = false;
      SetWindow ( 1 );
      redraw ();
      SetWindow ( 0 );
    }
    NLFinalSurfValid = swDisplayNLFinalPatches = false;
    if ( swDisplayFinalPatches )
      TurnFinalSurface ();
    break;

  default:
    break;
    }
    break;

case 1:
    switch ( er->id ) {
  case 3:  case 4:  case 5:
      RedrawGeomWindows ();
      er->redraw ( er );
      break;

  case 6:
      TurnFinalSurface ();
      er->redraw ( er );
      break;

  case 11:
      TurnNLFinalSurface ();
      er->redraw ( er );
      break;

  case 25:
      TurnConstraints ();
      swDisplayNLFinalPatches = false;
      PictureIsOn = false;
      er->redraw ( er );
      RedrawGeomWindows ();
      break;

  case 26:
      er->redraw ( er );
      break;

  case 27:
      if ( swZeroDer )
        menurect1b[1].redraw = DrawEmpty;
      else
        menurect1b[1].redraw = DrawSwitch;
      edrect1[6].redraw ( &edrect1[6] );
      break;

  case 28:
      er->redraw ( er );
      if ( swConstraintsOn ) {
        FinalSurfValid = false;
        RedrawGeomWindows ();
      }
      break;

  case 31:
      if ( rswGaussian )
        rswMean = false;
      edrect1[6].redraw ( &edrect1[6] );
      break;

  case 32:
      if ( rswMean )
        rswGaussian = false;
      edrect1[6].redraw ( &edrect1[6] );
      break;

  case 34:
      if ( rswSections )
        rswVDepRefl = rswVIndRefl = false;
      edrect1[6].redraw ( &edrect1[6] );
      break;

  case 35:
      if ( rswVDepRefl )
        rswSections = rswVIndRefl = false;
      edrect1[6].redraw ( &edrect1[6] );
      break;

  case 36:
      if ( rswVIndRefl )
        rswSections = rswVDepRefl = false;
      edrect1[6].redraw ( &edrect1[6] );
      break;

  case 37:
      edrect1[6].redraw ( &edrect1[6] );
      break;

  case 41:
      if ( eswSections )
        eswLines1 = eswLines2 = eswLight = false;
      eswEdRendering = eswSections || eswLines1 || eswLines2 || eswLight;
      redraw ();
      break;

  case 42:
      if ( eswLines1 )
        eswSections = eswLines2 = eswLight = false;
      eswEdRendering = eswSections || eswLines1 || eswLines2 || eswLight;
      redraw ();
      break;

  case 43:
      if ( eswLines2 )
        eswSections = eswLines1 = eswLight = false;
      eswEdRendering = eswSections || eswLines1 || eswLines2 || eswLight;
      redraw ();
      break;

  case 44:
      if ( eswLight )
        eswSections = eswLines1 = eswLines2 = false;
      eswEdRendering = eswSections || eswLines1 || eswLines2 || eswLight;
      redraw ();
      break;

  case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107:
  case 108: case 109: case 110: case 111: case 112: case 113: case 114: case 115:
  case 116: case 117: case 118: case 119: case 120: case 121: case 122: case 123:
  case 124: case 125: case 126: case 127: case 128: case 129: case 130: case 131:
  case 132:
      SwitchAConstraint ( er->id-100 );
      swDisplayNLFinalPatches = false;
      if ( swConstraintsOn ) {
        PictureIsOn = false;
        redraw ();
      }
      else
        edrect1[6].redraw ( &edrect1[6] );  /* more than one switch may change, */
                                            /* entire menu is to be updated */
      break;

  default: break;
    }
    break;

case 2:
    break;

default:
    break;
  }
} /*MySwitchProc*/

void GetSwitchData ( ed_rect *er, char **title, boolean **switchvar,
                     void (**switchproc)(ed_rect *er) )
{
  *switchproc = &MySwitchProc;
  *title = "";
  switch ( CurrentWindow () ) {
case 0:
    switch ( er->id ) {
  case  1: *switchvar = &HoleKSwitch[0];  break;
  case  2: *switchvar = &HoleKSwitch[1];  break;
  case  3: *switchvar = &HoleKSwitch[2];  break;
  case  4: *switchvar = &HoleKSwitch[3];  break;
  case  8: *switchvar = &swDisplayFirstPartition;  *title = "domain lines";    break;
  case  9: *switchvar = &swDisplayDomainCP;        *title = "control net";     break;
  case 10: *switchvar = &swDisplayDomNumbers;      *title = "numbers";         break;
  case 11: *switchvar = &swDisplayDomSurrPatches;  *title = "surr. patches";   break;
  case 12: *switchvar = &swDisplayDomCurves;       *title = "domain curves";   break;
  case 13: *switchvar = &swDisplayDomAuxPatches;   *title = "aux. patches";    break;
  case 14: *switchvar = &swDisplayDomPatches;      *title = "domain patches";  break;
  case 18: *switchvar = &swDisplayCentralPoint;    *title = "central point";   break;
  case 19: *switchvar = &swUseDerivatives1;        *title = "derivatives 1";   break;
  case 21: *switchvar = &swAltDomCurves;           *title = "alt. curves";     break;
  case 22: *switchvar = &swRestrictBasis;          *title = "restrict space";  break;
  case 23: *switchvar = &swCoonsPatches;           *title = "Coons space";     break;
  case 24: *switchvar = &swBezierPatches;          *title = "Bezier space";    break;
  default: break;
    }
    break;

case 1:
    switch ( er->id ) {
  case  3: *switchvar = &swDisplaySurfCP;          *title = "control net";     break;
  case  4: *switchvar = &swDisplaySurfNumbers;     *title = "numbers";         break;
  case  5: *switchvar = &swDisplaySurfPatches;     *title = "surface";         break;
  case  6: *switchvar = &swDisplayFinalPatches;    *title = "final patches A"; break;
  case 11: *switchvar = &swDisplayNLFinalPatches;  *title = "final patches B"; break;

  case 25: *switchvar = &swConstraintsOn;          *title = "constraints on";  break;
  case 26: *switchvar = &swPointMode;              *title = "point mode";      break;
  case 27: *switchvar = &swZeroDer;                *title = "zero";            break;
  case 28: *switchvar = &swNormalConstr;           *title = "normal";          break;
  case 100:
       if ( swZeroDer )
         *switchvar = &swZConstraint[0];
       else
         *switchvar = &swConstraint[0];
       *title = "central point";
       break;
  case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108:
  case 109: case 110: case 111: case 112: case 113: case 114: case 115: case 116:
  case 117: case 118: case 119: case 120: case 121: case 122: case 123: case 124:
  case 125: case 126: case 127: case 128: case 129: case 130: case 131: case 132:
      if ( swZeroDer )
        *switchvar = &swZConstraint[er->id-100];
      else
        *switchvar = &swConstraint[er->id-100];
      *title = "";
      break;

  case 31: *switchvar = &rswGaussian;  *title = "Gaussian";    break;
  case 32: *switchvar = &rswMean;      *title = "mean";        break;
  case 34: *switchvar = &rswSections;  *title = "sections";    break;
  case 35: *switchvar = &rswVDepRefl;  *title = "v.dep.refl";  break;
  case 36: *switchvar = &rswVIndRefl;  *title = "v.ind.refl";  break;
  case 37: *switchvar = &rswChess;     *title = "chess";       break;
  case 41: *switchvar = &eswSections;  *title = "sections";    break;
  case 42: *switchvar = &eswLines1;    *title = "v.dep.refl";  break;
  case 43: *switchvar = &eswLines2;    *title = "v.ind.refl";  break;
  case 44: *switchvar = &eswLight;     *title = "light";       break;
  default: break;
    }
    break;

case 2:
    break;

default:
    break;
  }
} /*GetSwitchData*/

void MyButtonProc ( ed_rect *er )
{
  switch ( CurrentWindow () ) {
case 0:
    switch ( er->id ) {
 case 0 :
     done = 1;
     break;

 case 15:
      PictureIsOn = false;
      FitToSurface ();
      break;

 case 16: case 20:
      SwitchMenu0 ();
      break;

 case 17:
      FixDomainCP ();
      redraw ();
      break;

 case 21:
      WriteInfo ();
      break;

 case 22:
      WriteFile ();
      break;

 default:
      break;
    }
    break;

case 1:
    switch ( er->id ) {
 case  8:
      Flatten ();
      break;

 case  9:
      FitToDomain ();
      break;

 case 10:
      SwitchMenu1 ( 1 );
      break;

 case 12:
      SwitchMenu1 ( 2 );
      break;

 case 38:
      SwitchMenu1 ( 0 );
      break;

 case 13:
      SetDefaultConstraints ();
      RedrawGeomWindows ();
      break;

 case 39:
      OnOffRendering ();
      break;

 default:
      break;
    }
    break;

case 2:
    break;

default:
    break;
  }
} /*MyButtonProc*/

void GetButtonData ( ed_rect *er, char **title,
                     void (**buttonproc)(ed_rect *er) )
{
  *buttonproc = &MyButtonProc;
  switch ( CurrentWindow () ) {
case 0:
    switch ( er->id ) {
 case  0: *title = "Quit";            break;
 case 15: *title = "Fit to surface";  break;
 case 16: *title = "Advanced";        break;
 case 17: *title = "Fix domain";      break;
 case 20: *title = "Basic";           break;
 case 21: *title = "Information";     break;
 case 22: *title = "Write";           break;
 default: break;
    }
    break;

case 1:
    switch ( er->id ) {
 case  8: *title = "Flatten";        break;
 case  9: *title = "Fit to domain";  break;
 case 10: *title = "Constraints";    break;
 case 12: *title = "Picture";        break;
 case 38: *title = "Data";           break;
 case 39:
      if ( RenderingIsOn )
        *title = "Stop";
      else
        *title = "Render";
      break;
 case 13: *title = "Current";        break;
 default: break;
    }
    break;

case 2:
    break;

default:
    break;
  }
} /*GetButtonData*/

void MySlideProc ( ed_rect *er )
{
  switch ( CurrentWindow () ) {
case 0:
    switch ( er->id ) {
  case 6:
  case 7:
  case 8:  SetDomainParam ();  break;
  default: break;
    }
    break;

case 1:
    switch ( er->id ) {
  case 0:
  case 1:
  case 2:
  case 3:  SetSurfParam ();  break;
  default: break;
    }
    break;

case 2:
    break;

default:
    break;
  }
} /*MySlideProc*/

void GetSlidebarData ( ed_rect *er, float **slipos,
                       void (**slideproc)(ed_rect *er) )
{
  *slideproc = &MySlideProc;
  switch ( CurrentWindow () ) {
case 0:
    switch ( er->id ) {
  case 6: *slipos = &domparams[0];  break;
  case 7: *slipos = &domparams[1];  break;
  case 8: *slipos = &domparams[2];  break;
  default: break;
    }
    break;

case 1:
    switch ( er->id ) {
  case 0: *slipos = &surfparams[0];  break;
  case 1: *slipos = &surfparams[1];  break;
  case 2: *slipos = &surfparams[2];  break;
  case 3: *slipos = &surfparams[3];  break;
  default: break;
    }
    break;

case 2:
    break;
  }
} /*GetSlidebarData*/

void GetMenuData ( ed_rect *er, int *menu_rect_num, ed_rect **menu )
{
  switch ( CurrentWindow () ) {
case 0:
    switch ( menu0 ) {
  case 0:
      *menu_rect_num = MENU0A_RECT_NUM;
      *menu = &menurect0a[0];
      break;
  case 1:
      *menu_rect_num = MENU0B_RECT_NUM;
      *menu = &menurect0b[0];
      break;
    }
    break;

case 1:
    switch ( menu1 ) {
  case 0:
      *menu_rect_num = MENU1A_RECT_NUM;
      *menu = &menurect1a[0];
      break;
  case 1:
      *menu_rect_num = MENU1B_RECT_NUM;
      *menu = &menurect1b[0];
      break;
  case 2:
      *menu_rect_num = MENU1C_RECT_NUM;
      *menu = &menurect1c[0];
      break;
    }
    break;

default:
    break;
  }
} /*GetMenuData*/

static char *InfoMsg[8] =
  { "Filling a polygonal hole in a piecewise bicubic surface with",
    "tangent plane continuity, by solving the biharmonic equation.",
    "",
    "This program is obsolete, try using the program policz.",
    "",
    "This program is a part of the BSTools package.",
    "(C) Copyright by Przemyslaw Kiciak, 2005,2007",
    NULL };

void ShowInfoScreen ()
{
  DisplayInfoMessage ( InfoMsg );
} /*ShowInfoScreen*/

void init_edwin ( int argc, char *argv[] )
{
  pkv_InitScratchMem ( SCRATCH_MEM_SIZE );

  SetupEdRect ( edrect0, 0, 0, WIDTH-110, HEIGHT-51, 110, 0,
                DomWinMsg, DrawDomWin );
  SetupEdRect ( edrect0, 1, 1, WIDTH-110, 50, 110, HEIGHT-50,
                KnotWinMsg, DrawKnotWin );
  SetupEdRect ( edrect0, 2, 2, 110, HEIGHT, 0, 0, MenuMsg, DrawMenu );
  SetupEdRect ( edrect0, 3, 3, 0, 0, 0, 0, ThinkingMsg, DrawEmpty );

  SetupEdRect ( menurect0a,  0,  0, 90, 19, 0, HEIGHT-19, ButtonMsg, DrawButton );
  SetupEdRect ( menurect0a,  1,  1, 16, 16,  0, 25, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0a,  2,  2, 16, 16, 20, 25, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0a,  3,  3, 16, 16, 40, 25, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0a,  4,  4, 16, 16, 60, 25, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0a,  5,  5, 21, 18, 80, 24, EmptyMsg, DisplayHoleK );
  SetupEdRect ( menurect0a,  6,  6, 109, 10, 0, 48, SlidebarMsg, DrawSlidebar );
  SetupEdRect ( menurect0a,  7,  7, 109, 10, 0, 66, SlidebarMsg, DrawSlidebar );
  SetupEdRect ( menurect0a,  8,  8, 109, 10, 0, 84, SlidebarMsg, DrawSlidebar );
  SetupEdRect ( menurect0a,  9,  8, 16, 16, 0, 148, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0a, 10,  9, 16, 16, 0, 168, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0a, 11, 10, 16, 16, 0, 188, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0a, 12, 11, 16, 16, 0, 208, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0a, 13, 12, 16, 16, 0, 228, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0a, 14, 13, 16, 16, 0, 248, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0a, 15, 14, 16, 16, 0, 268, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0a, 16, 17, 90, 19, 0, 101, ButtonMsg, DrawButton );
  SetupEdRect ( menurect0a, 17, 15, 90, 19, 0, 122, ButtonMsg, DrawButton );
  SetupEdRect ( menurect0a, 18, 16, 90, 19, 0,   0, ButtonMsg, DrawButton );
  SetupEdRect ( menurect0a, 19, 21, 90, 19, 0, HEIGHT-61, ButtonMsg, DrawButton );
  SetupEdRect ( menurect0a, 20, 22, 90, 19, 0, HEIGHT-40, ButtonMsg, DrawButton );

  SetupEdRect ( menurect0b,  0, 20, 90, 19, 0,         0, ButtonMsg, DrawButton );
  SetupEdRect ( menurect0b,  1, 18, 16, 16, 0,  25, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0b,  2, 19, 16, 16, 0,  45, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0b,  3, 21, 16, 16, 0,  65, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0b,  4, 22, 16, 16, 0,  85, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0b,  5, 23, 16, 16, 0, 105, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect0b,  6, 24, 16, 16, 0, 125, SwitchMsg, DrawSwitch );

  SetWinEdRect ( MAIN0_RECT_NUM, edrect0 );

  NewWindow ( "" );
  SetWindow ( 1 );

  SetupEdRect ( edrect1,  0, 0, 0, 0, 0, 0, EmptyMsg, DrawEmpty );
  SetupEdRect ( edrect1,  1, 1, 0, 0, 0, 0, ParWindowMsg, DrawParWindow );
  SetupEdRect ( edrect1,  2, 2, 0, 0, 0, 0, ParWindowMsg, DrawParWindow );
  SetupEdRect ( edrect1,  3, 3, 0, 0, 0, 0, ParWindowMsg, DrawParWindow );
  SetupEdRect ( edrect1,  4, 4, 0, 0, 0, 0, PerspWindowMsg, DrawPerspWindow );
  SetupEdRect ( edrect1,  5, 5, 0, 0, 0, 0, SpecialMsg, DrawSpecialWindow );
  SetupEdRect ( edrect1,  6, 6, 110, HEIGHT, 0, 0, MenuMsg, DrawMenu );

  SetupEdRect ( menurect1a, 0, 0, 109, 10, 0, 25, SlidebarMsg, DrawSlidebar );
  SetupEdRect ( menurect1a, 1, 1, 109, 10, 0, 43, SlidebarMsg, DrawSlidebar );
  SetupEdRect ( menurect1a, 2, 2, 109, 10, 0, 61, SlidebarMsg, DrawSlidebar );
  SetupEdRect ( menurect1a, 3, 3, 109, 10, 0, 79, EmptyMsg, DrawEmpty ); /* a slidebar only for a hexagonal hole */
  SetupEdRect ( menurect1a, 4, 3, 16, 16, 0, 142, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1a, 5, 4, 16, 16, 0, 162, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1a, 6, 5, 16, 16, 0, 182, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1a, 7, 6, 16, 16, 0, 202, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1a, 8, 8, 90, 19, 0, 97, ButtonMsg, DrawButton );
  SetupEdRect ( menurect1a, 9, 9, 90, 19, 0, 118, ButtonMsg, DrawButton );
  SetupEdRect ( menurect1a, 10, 10, 90, 19, 0, 0, ButtonMsg, DrawButton );
  SetupEdRect ( menurect1a, 11, 11, 16, 16, 0, 222, SwitchMsg, DrawSwitch );

  SetupEdRect ( menurect1b,  0, 12, 90, 19, 0, 0, ButtonMsg, DrawButton );
  SetupEdRect ( menurect1b,  1, 100, 16, 16,  0, 125, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b,  2, 101, 16, 16, 20, 145, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b,  3, 102, 16, 16, 40, 145, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b,  4, 103, 16, 16, 60, 145, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b,  5, 104, 16, 16, 80, 145, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b,  6, 105, 16, 16, 20, 165, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b,  7, 106, 16, 16, 40, 165, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b,  8, 107, 16, 16, 60, 165, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b,  9, 108, 16, 16, 80, 165, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 10, 109, 16, 16, 20, 185, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 11, 110, 16, 16, 40, 185, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 12, 111, 16, 16, 60, 185, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 13, 112, 16, 16, 80, 185, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 14, 113, 16, 16, 20, 205, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 15, 114, 16, 16, 40, 205, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 16, 115, 16, 16, 60, 205, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 17, 116, 16, 16, 80, 205, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 18, 117, 16, 16, 20, 225, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 19, 118, 16, 16, 40, 225, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 20, 119, 16, 16, 60, 225, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 21, 120, 16, 16, 80, 225, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 22, 121, 16, 16, 20, 245, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 23, 122, 16, 16, 40, 245, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 24, 123, 16, 16, 60, 245, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 25, 124, 16, 16, 80, 245, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 26, 125, 16, 16, 20, 265, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 27, 126, 16, 16, 40, 265, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 28, 127, 16, 16, 60, 265, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 29, 128, 16, 16, 80, 265, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 30, 129, 16, 16, 20, 285, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 31, 130, 16, 16, 40, 285, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 32, 131, 16, 16, 60, 285, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 33, 132, 16, 16, 80, 285, EmptyMsg, DrawEmpty );
  SetupEdRect ( menurect1b, 34,  13, 90, 19,  0,  21, ButtonMsg, DrawButton );
  SetupEdRect ( menurect1b, 35,  25, 16, 16,  0,  45, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 36,  26, 16, 16,  0,  65, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 37,  27, 16, 16,  0, 105, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1b, 38,  28, 16, 16,  0,  85, SwitchMsg, DrawSwitch );

  SetupEdRect ( menurect1c,  0,  30, 109, 16, 0, 20, EmptyMsg, DrawText );
  SetupEdRect ( menurect1c,  1,  31, 16, 16, 0,  36, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1c,  2,  32, 16, 16, 0,  56, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1c,  3,  33, 16, 16, 0,  76, EmptyMsg, DrawText );   
  SetupEdRect ( menurect1c,  4,  34, 16, 16, 0,  96, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1c,  5,  35, 16, 16, 0, 116, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1c,  6,  36, 16, 16, 0, 136, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1c,  7,  37, 16, 16, 0, 156, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1c,  8,  38, 90, 18, 0,   0, ButtonMsg, DrawButton );
  SetupEdRect ( menurect1c,  9,  39, 90, 18, 0, HEIGHT-19, ButtonMsg, DrawButton );
  SetupEdRect ( menurect1c, 10,  40, 109, 16, 0, 176, EmptyMsg, DrawText );
  SetupEdRect ( menurect1c, 11,  41, 16, 16, 0, 196, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1c, 12,  42, 16, 16, 0, 216, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1c, 13,  43, 16, 16, 0, 236, SwitchMsg, DrawSwitch );
  SetupEdRect ( menurect1c, 14,  44, 16, 16, 0, 256, SwitchMsg, DrawSwitch );

  SetWinEdRect ( MAIN1_RECT_NUM, edrect1 );


  SetupRefBox ( -1.0, 1.0, -1.0, 1.0, -1.0, 1.0 );
  InitProjections ();
  InitSurface ( 5 );
  SetConstraintSwitches ( 5 );
  ResizeWindow0 ();
  ResizeWindow1 ( true );
  InitRenderer ();

  redraw_all ();
  SetWindow ( 0 );
  ShowInfoScreen ();
} /*init_edwin*/

void destroy_edwin ()
{
  printf ( "Scratch memory used: %d out of %d bytes\n",
           (int)pkv_MaxScratchTaken(), SCRATCH_MEM_SIZE );
  DestroyRenderer ();
  pkv_DestroyScratchMem ();
} /*destroy_edwin*/

void resize_edwin ()
{
  GetWindowSize ();
  switch ( CurrentWindow () ) {
case 0: ResizeWindow0 ();  break;
case 1: ResizeWindow1 ( false );  break;
case 2: break;
default: break;
  }
  redraw ();
} /*resize_edwin*/

void process_key ( unsigned int *key )
{
  int dx, dy;

  switch ( *key ) {
case 'M':
    XResizeWindow ( thedisplay, thewindow, MAX_WIDTH, MAX_HEIGHT );
    break;

case 'm':
    XResizeWindow ( thedisplay, thewindow, WIDTH, HEIGHT );
    break;

case 'D':
case 'd':
    XMoveWindow ( thedisplay, thewindow, 0, 512 );
    break;

case 'Q':
case 'q':
    done = 1;
    break;

default:
    switch ( thekeysym ) {
  case 0xFF51: /* key left */
      dx = -1;  dy = 0;
      goto cont;
  case 0xFF52: /* key up */
      dx = 0;  dy = -1;
      goto cont;
  case 0xFF53: /* key right */
      dx = +1;  dy = 0;
      goto cont;
  case 0xFF54: /* key down */
      dx = 0;  dy = +1;
  cont:
      XWarpPointer ( thedisplay, None, None, 0, 0, 0, 0, dx, dy );
      dispatch_message ( msg_MMOVE, mouse_buttons, mouse_x+dx, mouse_y+dy );
      break;
  default:
      break;
    }
    return;
  }
  *key = 0;
} /*process_key*/

