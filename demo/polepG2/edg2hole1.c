
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
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
#include "eg2holef.h"

#include "oldxgedit.h"
#include "datagenf.h"
#include "g2ekernel.h"
#include "edg2hole.h"
#include "render.h"

int hole_k = 5;
boolean HoleKSwitch[4] = {false,true,false,false};

boolean swDisplayFirstPartition = false;
boolean swDisplayDomainCP = true;
boolean swDisplayDomSurrPatches = true;
boolean swDisplayDomCurves = false;
boolean swDisplayDomAuxPatches = false;
boolean swDisplayDomPatches = false;
boolean swDisplayPartition = false;
boolean swDisplayDomNumbers = false;
boolean swDisplayCentralPoint = false;
boolean swUseDerivatives1 = false;
boolean swAltDomCurves = false;
boolean swRestrictBasis = false;

static float slDomainWinDiag;
static int current_np;
static int ckn_i, ckn_j;

/* ///////////////////////////////////////////////////////////////////////// */
void SetHoleK ( int k )
{
  SetConstraintSwitches ( k );
  swDisplayCentralPoint = swUseDerivatives1 = swAltDomCurves = false;
  HoleKSwitch[0] = HoleKSwitch[1] = HoleKSwitch[2] = HoleKSwitch[3] = false;
  switch ( k ) {
case 3: HoleKSwitch[0] = true;  break;
case 5: HoleKSwitch[1] = true;  break;
case 6: HoleKSwitch[2] = true;  break;
case 8: HoleKSwitch[3] = true;  break;
default: exit ( 1 );
  }
  InitSurface ( hole_k = k );
  ProjectDomainCP ();
  ProjectSurfCP ();
  RecreateDomain ();
  swDisplayCentralPoint   = false;
  swDisplayNLFinalPatches = false;
  swDisplayFinalPatches   = FinalSurfValid   = false;
  swConstraintsOn         = ConstraintsValid = false;
  ConstraintsValid        = false;
  ResetRendLines ();
} /*SetHoleK*/

void SetDomZoom ( float zf )
{
  float f;

  f = slDomainWinDiag*exp ( zf );
  f = max ( f, MIN_DOMZOOM );
  f = min ( f, MAX_DOMZOOM );
  DomPPos.vd.para.diag = f;
  CameraSetMappingf ( &DomPPos );
} /*SetDomZoom*/

void DomWinMsg ( ed_rect *er, int msg, int key, int x, int y )
{
  switch ( stan ) {
case STATE_NOTHING:
    switch ( msg ) {
  case msg_MCLICK:
      if ( (key & mouse_LBUTTON_DOWN) && (key & mouse_LBUTTON_CHANGE) ) {
        if ( swDisplayCentralPoint ) {
          current_np = FindNearestCPointPart ( x, y );
          if ( current_np >= 0 ) {
            xx = x;  yy = y;
            stan = STATE_MOVING_CENTP;
            focus = er;
          }
        }
        else if ( swDisplayFirstPartition ) {
          current_np = FindNearestFirstPartitionVector ( x, y );
          if ( current_np >= 0 ) {
            xx = x;  yy = y;
            stan = STATE_MOVING_FIRSTP;
            focus = er;
          }
        }
        else if ( swDisplayDomainCP ) {
          current_np = FindNearestDomCP ( x, y );
          if ( current_np >= 0 ) {
            xx = x;  yy = y;
            stan = STATE_MOVING_DOMCP;
            focus = er;
          }
        }
      }
      else if ( (key & mouse_RBUTTON_DOWN) && (key & mouse_RBUTTON_CHANGE) ) {
        slDomainWinDiag = DomPPos.vd.para.diag;
        xx = yy = y;
        stan = STATE_DOMZOOM;
        focus = er;
      }
    }
    break;

case STATE_MOVING_DOMCP:
    if ( msg == msg_MMOVE || msg == msg_MCLICK ) {
      if ( !(key & mouse_LBUTTON_DOWN) ) {
        stan = STATE_NOTHING;
        focus = NULL;
      }
      else {
        BoundPoint ( er, &x, &y );
        if ( x != xx || y != yy ) {
          SetDomCP ( current_np, x, y );
          er->redraw ( er );
          PictureIsOn = false;
          if ( RenderingIsOn )
            OffRendering ();
          if ( swDisplayFinalPatches ) {
            swDisplayFinalPatches = false;
            SetWindow ( 1 );
            redraw ();
            SetWindow ( 0 );
          }
          xx = x; yy = y;
        }
      }
    }
    break;

case STATE_MOVING_FIRSTP:
    if ( msg == msg_MMOVE || msg == msg_MCLICK ) {
      if ( !(key & mouse_LBUTTON_DOWN) ) {
        stan = STATE_NOTHING;
        focus = NULL;
      }
      else {
        BoundPoint ( er, &x, &y );
        if ( x != xx || y != yy ) {
          if ( SetFirstPartitionVector ( current_np, x, y ) ) {
            SetDomainParam ();
            er->redraw ( er );
            xx = x; yy = y;
          }
        }
      }
    }
    break;

case STATE_MOVING_CENTP:
    if ( msg == msg_MMOVE || msg == msg_MCLICK ) {
      if ( !(key & mouse_LBUTTON_DOWN) ) {
        stan = STATE_NOTHING;
        focus = NULL;
      }
      else {
        BoundPoint ( er, &x, &y );
        if ( x != xx || y != yy ) {
          if ( SetCPointPart ( current_np, x, y ) ) {
            RecreateDomain ();
            if ( swDisplayFinalPatches ) {
              swDisplayFinalPatches = false;
              SetWindow ( 1 );
              redraw ();
              SetWindow ( 0 );
            }
            er->redraw ( er );
            xx = x;  yy = y;
          }
        }
      }
    }
    break;

case STATE_DOMZOOM:
    switch ( msg ) {
  case msg_MMOVE:
      if ( !(key & mouse_RBUTTON_DOWN) ) {
        stan = STATE_NOTHING;
        focus = NULL;
      }
      else if ( y != yy ) {
        SetDomZoom ( (float)(y-xx)/(float)er->h );
        ProjectDomainCP ();
        er->redraw ( er );
        yy = y;
      }
      break;

  case msg_MCLICK:
      if ( !(key & mouse_RBUTTON_DOWN) ) {
        stan = STATE_NOTHING;
        focus = NULL;
      }
      break;

  default:
      break;
    }
    break;

default:
    break;
  }
} /*DomWinMsg*/

void DrawDomWin ( ed_rect *er )
{
  XSetForeground ( thedisplay, thegc, c_dk_blue );
  XFillRectangle ( thedisplay, thepixmap, thegc,
                   er->x+1, er->y+1, er->w-2, er->h-2 );
  XSetForeground ( thedisplay, thegc, c_lt_blue );
  XDrawRectangle ( thedisplay, thepixmap, thegc,
                   er->x, er->y, er->w-1, er->h-1 );

  if ( swDisplayDomSurrPatches )
    RedrawDomainSurrPatches ();
  if ( swDisplayDomCurves & !swDisplayDomPatches )
    RedrawDomainCurves ();
  if ( swDisplayPartition )
    RedrawDomainPartition ();
  if ( swDisplayDomAuxPatches )
    RedrawDomainAuxPatches ();
  if ( swDisplayDomPatches )
    RedrawDomainPatches ();
  if ( swDisplayDomainCP )
    RedrawDomainCP ( !swDisplayFirstPartition && !swDisplayCentralPoint );
  if ( swDisplayFirstPartition )
    RedrawFirstPartition ( !swDisplayCentralPoint );
  if ( swDisplayDomNumbers )
    RedrawDomainNumbers ();
  if ( swDisplayCentralPoint )
    RedrawDomainCentralPoint ();

  XCopyArea ( thedisplay, thepixmap, thewindow, thegc,
              er->x, er->y, er->w, er->h, er->x, er->y );
} /*DrawDomWin*/

void KnotWinMsg ( ed_rect *er, int msg, int key, int x, int y )
{
  switch ( stan ) {
case STATE_NOTHING:
    if ( msg == msg_MCLICK ) {
      if ( (key & mouse_LBUTTON_DOWN) && (key & mouse_LBUTTON_CHANGE) ) {
        if ( FindNearestKnot ( er->w, er->h, er->x, er->y,
                               x, y, &ckn_i , &ckn_j ) ) {
          xx = x;
          stan = STATE_MOVINGKNOT;
          focus = er;
        }
      }
    }
    break;

case STATE_MOVINGKNOT:
    if ( msg == msg_MMOVE ) {
      if ( key & mouse_LBUTTON_DOWN ) {
        BoundPoint ( er, &x, &y );
        if ( x != xx ) {
          if ( SetKnot ( er->w, er->h, er->x, er->y, ckn_i, ckn_j, x ) ) {
            xx = x;
            er->redraw ( er );
            edrect0[0].redraw ( &edrect0[0] );
            PictureIsOn = false;
            if ( RenderingIsOn )
              OffRendering ();
            if ( swDisplayFinalPatches ) {
              swDisplayFinalPatches = false;
              SetWindow ( 1 );
              redraw ();
              SetWindow ( 0 );
            }
            else
              RedrawGeomWindows ();
          }
        }
      }
      else {
        stan = STATE_NOTHING;
        focus = NULL;
      }
    }
    else if ( msg == msg_MCLICK ) {
      if ( !(key & mouse_LBUTTON_DOWN) ) {
        stan = STATE_NOTHING;
        focus = NULL;
      }
    }
    break;

default:
    break;
  }
} /*KnotWinMsg*/

void DrawKnotWin ( ed_rect *er )
{
  XSetForeground ( thedisplay, thegc, c_dk_blue );
  XFillRectangle ( thedisplay, thepixmap, thegc,
                   er->x+1, er->y+1, er->w-2, er->h-2 );
  XSetForeground ( thedisplay, thegc, c_lt_blue );
  XDrawRectangle ( thedisplay, thepixmap, thegc,
                   er->x, er->y, er->w-1, er->h-1 );

  RedrawHoleKnots ( er->w, er->h, er->x, er->y );

  XCopyArea ( thedisplay, thepixmap, thewindow, thegc,
              er->x, er->y, er->w, er->h, er->x, er->y );
} /*DrawKnotWin*/

void SetDomainParam ()
{
  int w;

  InitDomain ( hole_k, firstpartition, domparams, domcp );
  RecreateDomain ();
  ProjectDomainCP ();
  swDisplayFinalPatches   = FinalSurfValid   = false;
  swDisplayNLFinalPatches = NLFinalSurfValid = false;
  swConstraintsOn         = ConstraintsValid = false;
  edrect0[0].redraw ( &edrect0[0] );
  w = CurrentWindow ();
  PictureIsOn = false;
  if ( RenderingIsOn )
    OffRendering ();
  SetWindow ( 1 );
  redraw ();
  SetWindow ( w );
} /*SetDomainParam*/

void ResizeWindow0 ()
{
  int w, h;

  w = current_width-110;
  h = current_height-51;
  edrect0[0].w = w;  edrect0[0].h = h;
  edrect0[1].w = w;  edrect0[1].y = h+1;
  edrect0[2].h = current_height;
  menurect0a[0].y = current_height-19;
  menurect0a[19].y = current_height-61;
  menurect0a[20].y = current_height-40;

  SetupDomPPos ( w, h, edrect0[0].x, edrect0[0].y, DomPPos.vd.para.diag );
  ProjectDomainCP ();
} /*ResizeWindow0*/

void FitToSurface ()
{
  FitDomainToSurface ();
  swDisplayFinalPatches   = FinalSurfValid     = false;
  swDisplayNLFinalPatches = NLFinalSurfValid   = false;
  FormMatrixValid         = ExtFormMatrixValid = false;
  swConstraintsOn         = ConstraintsValid   = false;
  PictureIsOn = false;
  redraw_all ();
} /*FitToSurface*/

void FixDomainCP ()
{
  StoreDomain ( hole_k, domcp );
  domparams[2] = 1.0;
} /*FixDomainCP*/

void SwitchMenu0 ()
{
  menu0 = (menu0+1) % 2;
  edrect0[2].redraw ( &edrect0[2] );
} /*SwitchMenu0*/

