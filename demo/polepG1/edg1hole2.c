
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2007                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/times.h>
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

boolean swDisplaySurfCP         = true;
boolean swDisplaySurfPatches    = true;
boolean swDisplayFinalPatches   = false;
boolean swDisplaySurfNumbers    = false;
boolean swDisplayNLFinalPatches = false;

static float zzoom;
static float fwinxx = 0.5;
static float fwinyy = 0.5;
static int gwinxx, gwinyy;

static int current_point;

boolean PictureIsOn = false;
static clock_t lastupdate, ticks_per_second;

/* ////////////////////////////////////////////////////////////////////////// */
static int cursorwin;

void RedrawGeomWindows ()
{
  int i, w;

  w = CurrentWindow ();
  SetWindow ( 1 );
  for ( i = FIRST_GEOMW; i < FIRST_GEOMW+4; i++ )
    edrect1[i].redraw ( &edrect1[i] );
  SetWindow ( w );
} /*RedrawGeomWindows*/

void SetSurfParam ()
{
  hole_np = InitHole ( hole_k, surfparams, surfcp );
  FinalSurfValid = NLFinalSurfValid = PictureIsOn = false;
  swDisplayNLFinalPatches = false;
  if ( RenderingIsOn )
    OffRendering ();
  ProjectSurfCP ();
  redraw ();
} /*SetSurfParam*/

void TurnFinalSurface ()
{
  PictureIsOn = false;
  if ( RenderingIsOn )
    OffRendering ();
  if ( swDisplayFinalPatches ) {
    if ( !FinalSurfValid ) {
      if ( swCoonsPatches ) {
        if ( !FormMatrixValid ) {
          cursorwin = CurrentWindow ();
          XDefineCursor ( thedisplay, thewindow, thecursor1 );
          stan = STATE_THINKING0;
          focus = &edrect0[3];
          return;
        }
        if ( FormMatrixValid )
          UpdateFinalSurf ();
      }
      else {
        if ( !ExtFormMatrixValid ) {
          cursorwin = CurrentWindow ();
          XDefineCursor ( thedisplay, thewindow, thecursor1 );
          stan = STATE_THINKING0;
          focus = &edrect0[3];
          return;
        }
        if ( ExtFormMatrixValid )
          UpdateFinalSurf ();
      }
    }
    if ( !FinalSurfValid )
      swDisplayFinalPatches = false;
  }
  RedrawGeomWindows ();
} /*TurnFinalSurface*/

void TurnNLFinalSurface ()
{
  PictureIsOn = false;
  if ( RenderingIsOn )
    OffRendering ();
  if ( swDisplayNLFinalPatches && !NLFinalSurfValid ) {
    cursorwin = CurrentWindow ();
    XDefineCursor ( thedisplay, thewindow, thecursor1 );
    stan = STATE_THINKING1;
    focus = &edrect0[3];
    return;
  }
  RedrawGeomWindows ();
} /*TurnNLFinalSurface*/

void Flatten ()
{
  FlattenSurface ();
  FinalSurfValid = NLFinalSurfValid = PictureIsOn = false;
  if ( RenderingIsOn )
    OffRendering ();
  ProjectSurfCP ();
  TurnFinalSurface ();
} /*Flatten*/

void FitToDomain ()
{
  FitSurfaceToDomain ();
  FinalSurfValid = NLFinalSurfValid = PictureIsOn = false;
  if ( RenderingIsOn )
    OffRendering ();
  ProjectSurfCP ();
  TurnFinalSurface ();
} /*FitToDomain*/

/* ////////////////////////////////////////////////////////////////////////// */
void DrawParWindow ( ed_rect *er )
{
  int id;

  id = er->id;
  XSetForeground ( thedisplay, thegc, c_dk_blue );
  XFillRectangle ( thedisplay, thepixmap, thegc,
             er->x+1, er->y+1, er->w-2, er->h-2 );
  XSetForeground ( thedisplay, thegc, c_lt_blue );
  XDrawRectangle ( thedisplay, thepixmap, thegc,
             er->x, er->y, er->w-1, er->h-1 );

  if ( swDisplaySurfPatches )
    RedrawSurface ( id-FIRST_GEOMW );
  if ( swDisplayFinalPatches )
    RedrawFinalPatches ( id-FIRST_GEOMW );
  if ( swDisplayNLFinalPatches )
    RedrawNLFinalPatches ( id-FIRST_GEOMW );
  if ( swDisplaySurfCP )
    RedrawSurfaceCP ( id-FIRST_GEOMW, (!swConstraintsOn) && (!eswEdRendering) );
  if ( swConstraintsOn )
    RedrawConstraintsHandle ( id-FIRST_GEOMW, !eswEdRendering );
  if ( swDisplaySurfNumbers )
    RedrawSurfNumbers ( id-FIRST_GEOMW );
  if ( eswEdRendering )
    RedrawRendPoints ( id-FIRST_GEOMW );

  XCopyArea ( thedisplay, thepixmap, thewindow, thegc,
              er->x, er->y, er->w, er->h, er->x, er->y );
} /*DrawParWindow*/

void SetParZoom ( float zf )
{
  float f;

  f = zzoom*exp(zf);
  f = max ( f, MIN_PARZOOM );
  f = min ( f, MAX_PARZOOM );
  SetupRefBox ( -f, f, -f, f, -f, f );
  SetupParProj ( &edrect1[FIRST_GEOMW] );
  ProjectSurfCP ();
} /*SetParZoom*/

void ParWindowMsg ( ed_rect *er, int msg, int key, int x, int y )
{
  int i, id;

  id = er->id;
  switch ( stan ) {
case STATE_NOTHING:
    switch ( msg ) {
  case msg_MCLICK:
      if ( (swDisplaySurfCP || swConstraintsOn || eswEdRendering ) &&
           (key & mouse_LBUTTON_DOWN) && (key & mouse_LBUTTON_CHANGE) ) {
        if ( (current_point = FindNearestPoint ( id-FIRST_GEOMW, x, y )) >= 0 ) {
          stan = STATE_MOVINGPOINT;
          focus = er;
          xx = x;
          yy = y;
        }
      }
      else if ( (key & mouse_RBUTTON_DOWN) && (key & mouse_RBUTTON_CHANGE) ) {
        stan = STATE_PARZOOM;
        focus = er;
        xx = yy = y;
        zzoom = 0.5*(RefBBox.x1-RefBBox.x0);
      }
      break;
    }
    break;

case STATE_PARZOOM:
    switch ( msg ) {
case msg_MMOVE:
      if ( !(key & mouse_RBUTTON_DOWN) ) {
        stan = STATE_NOTHING;
        focus = NULL;
      }
      else if ( y != yy ) {
        SetParZoom ( (float)(y-xx)/(float)er->h );
        for ( i = FIRST_GEOMW; i < FIRST_GEOMW+3; i++ ) {
          ProjectScene ( i-FIRST_GEOMW );
          edrect1[i].redraw ( &edrect1[i] );
        }
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

case STATE_MOVINGPOINT:
    switch ( msg ) {
case msg_MMOVE:
      if ( key & mouse_LBUTTON_DOWN ) {
        BoundPoint ( er, &x, &y );
        if ( x != xx || y != yy ) {
          SetCPoint ( id-FIRST_GEOMW, current_point, x, y );
          PictureIsOn = false;
          if ( RenderingIsOn )
            OffRendering ();
          RedrawGeomWindows ();
          xx = x;
          yy = y;
        }
      }
      else {
        stan = STATE_NOTHING;
        focus = NULL;
      }
      break;

case msg_MCLICK:
      if ( !(key & mouse_LBUTTON_DOWN) ) {
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
} /*ParWindowMsg*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawPerspWindow ( ed_rect *er )
{
  int id;

  id = er->id;
  if ( !PictureIsOn ) {
    XSetForeground ( thedisplay, thegc, c_dk_blue );
    XFillRectangle ( thedisplay, thepixmap, thegc,
               er->x+1, er->y+1, er->w-2, er->h-2 );
    XSetForeground ( thedisplay, thegc, c_lt_blue );
    XDrawRectangle ( thedisplay, thepixmap, thegc,
               er->x, er->y, er->w-1, er->h-1 );

    if ( swDisplaySurfPatches )
      RedrawSurface ( id-FIRST_GEOMW );
    if ( swDisplayFinalPatches )
      RedrawFinalPatches ( id-FIRST_GEOMW );
    if ( swDisplayNLFinalPatches )
      RedrawNLFinalPatches ( id-FIRST_GEOMW );
    if ( swDisplaySurfCP )
      RedrawSurfaceCP ( id-FIRST_GEOMW, (!swConstraintsOn) && (!eswEdRendering) );
    if ( swConstraintsOn )
      RedrawConstraintsHandle ( id-FIRST_GEOMW, !eswEdRendering );
    if ( swDisplaySurfNumbers )
      RedrawSurfNumbers ( id-FIRST_GEOMW );
    if ( eswEdRendering )
      RedrawRendPoints ( id-FIRST_GEOMW );
  }
  else {
    XPutImage ( thedisplay, thepixmap, thegc, theimage, 0, 0,
                er->x, er->y, er->w, er->h );
  }
  XCopyArea ( thedisplay, thepixmap, thewindow, thegc,
              er->x, er->y, er->w, er->h, er->x, er->y );
} /*DrawPerspWindow*/

void SetZoom ( float zf )
{
  float f;

  f = zzoom*exp(zf);
  f = max ( f, MIN_ZOOM );
  f = min ( f, MAX_ZOOM );
  CameraSetFf ( &CPos, f );
} /*SetZoom*/

void PerspWindowMsg ( ed_rect *er, int msg, int key, int x, int y )
{
  int      id;
  vector3f v;
  float    a, angle;

  id = er->id;
  switch ( stan ) {

case STATE_TURNINGVIEWER:
    if ( msg == msg_MMOVE || msg == msg_MCLICK ) {
      if ( key & mouse_LBUTTON_DOWN ) {
        if ( x != xx || y != yy ) {
          SetVector3f ( &v, (y-yy)/CPos.yscale, (xx-x)/CPos.xscale, 0.0 );
          a = 1000.0*sqrt ( ((float)v.x*v.x + v.y*v.y)/
                    ((float)CPos.width*CPos.width + CPos.height*CPos.height) );
          angle = -(2.0*a+1.0)*a;
          CameraRotVCf ( &CPos, &v, angle );
          ProjectScene ( id-FIRST_GEOMW );
          PictureIsOn = false;
          if ( RenderingIsOn )
            OffRendering ();
          er->redraw ( er );
          xx = x;
          yy = y;
        }
      }
      else {
        stan = STATE_NOTHING;
        focus = NULL;
      }
    }
    break;

case STATE_ZOOM:
    if ( msg == msg_MMOVE || msg == msg_MCLICK ) {
      if ( key & mouse_RBUTTON_DOWN ) {
        if ( y != yy ) {
          SetZoom ( (float)(xx-y)/(float)edrect1[id].h );
          ProjectScene ( id-FIRST_GEOMW );
          PictureIsOn = false;
          if ( RenderingIsOn )
            OffRendering ();
          er->redraw ( er );
          yy = y;
        }
      }
      else {
        stan = STATE_NOTHING;
        focus = NULL;
      }
    }
    break;

case STATE_NOTHING:
    if ( msg == msg_MCLICK ) {
      if ( key & mouse_LBUTTON_DOWN ) {
        stan = STATE_TURNINGVIEWER;
        focus = er;
        xx = x;
        yy = y;
      }
      else if ( key & mouse_RBUTTON_DOWN ) {
        stan = STATE_ZOOM;
        focus = er;
        xx = yy = y;
        zzoom = CPos.vd.persp.f;
      }
    }
    else if ( msg == msg_KEY ) {
      if ( key == 'r' || key == 'R' ) {
        ResetCPos ();
        ProjectScene ( id-FIRST_GEOMW );
        PictureIsOn = false;
        if ( RenderingIsOn )
          OffRendering ();
        er->redraw ( er );
      }
    }
    break;
  }
} /*PerspWindowMsg*/

void RestrictGWinXY ()
{
  gwinxx = max ( gwinxx, 130 );
  gwinxx = min ( gwinxx, current_width-20 );
  gwinyy = max ( gwinyy, 20 );
  gwinyy = min ( gwinyy, current_height-20 );
} /*RestrictGWinXY*/

void CompFWinXY ()
{
  fwinxx = (float)(gwinxx-110)/(float)(current_width-110);
  fwinyy = (float)gwinyy/(float)current_height;
} /*CompFWinXY*/

void DrawSpecialWindow ( ed_rect *er )
{
  XSetForeground ( thedisplay, thegc, c_black );
  XFillRectangle ( thedisplay, thewindow, thegc,
      edrect1[FIRST_GEOMW+0].x+edrect1[FIRST_GEOMW+0].w, 0,
      3, edrect1[FIRST_GEOMW+2].y+edrect1[FIRST_GEOMW+2].h );
  XFillRectangle ( thedisplay, thewindow, thegc,
      edrect1[FIRST_GEOMW+0].x, edrect1[FIRST_GEOMW+2].y-3,
      edrect1[FIRST_GEOMW+0].w+edrect1[FIRST_GEOMW+1].w+3, 3 );
} /*DrawSpecialWindow*/

void SpecialMsg ( ed_rect *er, int msg, int key, int x, int y )
{
  boolean cx, cy;

  switch ( stan ) {
case STATE_NOTHING:
    switch ( msg ) {
  case msg_ENTERING:
      XDefineCursor ( thedisplay, thewindow, thecursor0 );
      break;

  case msg_EXITING:
      SetWindow ( 1 );
      XDefineCursor ( thedisplay, thewindow, None );
      break;

  case msg_MCLICK:
      if ( key & mouse_LBUTTON_DOWN ) {
        cx = abs ( x - gwinxx ) <= 10;
        cy = abs ( y - gwinyy ) <= 10;
        if ( cx && cy ) stan = STATE_RESIZE_XY;
        else if ( cx )  stan = STATE_RESIZE_X;
        else if ( cy )  stan = STATE_RESIZE_Y;
        if ( stan != STATE_NOTHING ) {
          focus = er;
          xx = x;
          yy = y;
        }
      }
      break;

  case msg_KEY:
      if ( key == 'r' || key == 'R' ) {
        fwinxx = fwinyy = 0.5;
        ResizeWindow1 ( false );
        RedrawGeomWindows ();
        er->redraw ( er );
      }
      break;

  default:
      break;
    }

case STATE_RESIZE_X:
    if ( msg == msg_MCLICK ) {
      if ( !(key & mouse_LBUTTON_DOWN) )
        goto exit_resizing;
    }
    else if ( msg == msg_MMOVE ) {
      if ( !(key & mouse_LBUTTON_DOWN) )
        goto exit_resizing;
      if ( x != xx ) {
        gwinxx = xx = x;
        goto resize;
      }
    }
    break;

case STATE_RESIZE_Y:
    if ( msg == msg_MCLICK ) {
      if ( !(key & mouse_LBUTTON_DOWN) )
        goto exit_resizing;
    }
    else if ( msg == msg_MMOVE ) {
      if ( !(key & mouse_LBUTTON_DOWN) )
        goto exit_resizing;
      if ( y != yy ) {
        gwinyy = yy = y;
        goto resize;
      }
    }
    break;

case STATE_RESIZE_XY:
    if ( msg == msg_MCLICK ) {
      if ( !(key & mouse_LBUTTON_DOWN) )
        goto exit_resizing;
    }
    else if ( msg == msg_MMOVE ) {
      if ( !(key & mouse_LBUTTON_DOWN) ) {
exit_resizing:
        stan = STATE_NOTHING;
        focus = NULL;
        return;
      }
      if ( x != xx || y != yy ) {
        gwinxx = xx = x;
        gwinyy = yy = y;
resize:
        RestrictGWinXY ();
        CompFWinXY ();
        ResizeWindow1 ( false );
        RedrawGeomWindows ();
        er->redraw ( er );
      }
    }
    break;

default:
    ;
  }
} /*SpecialMsg*/

void ThinkingMsg ( ed_rect *er, int msg, int key, int x, int y )
{
  char *errmsg;

  switch ( stan ) {
case STATE_THINKING0:
    if ( swCoonsPatches ) {
      UpdateFormMatrix ();
      if ( FormMatrixValid )
        UpdateFinalSurf ();
    }
    else {
      UpdateExtFormMatrix ();
      if ( ExtFormMatrixValid )
        UpdateFinalSurf ();
    }
    SetWindow ( cursorwin );
    XDefineCursor ( thedisplay, thewindow, None );
    SetWindow ( 1 );
    if ( FinalSurfValid )
      redraw ();
    else {
      swDisplayFinalPatches = false;
      g1h_GetErrorCodef ( domain, &errmsg );
      DisplayErrorMessage ( errmsg );
      return;
    }
    break;

case STATE_THINKING1:
    UpdateNLFinalSurf ();
    SetWindow ( cursorwin );
    XDefineCursor ( thedisplay, thewindow, None );
    SetWindow ( 1 );
    if ( NLFinalSurfValid )
      redraw ();
    else {
      swDisplayNLFinalPatches = false;
      g1h_GetErrorCodef ( domain, &errmsg );
      DisplayErrorMessage ( errmsg );
      return;
    }
    break;

default:
    break;
  }
  stan = STATE_NOTHING;
  focus = NULL;
} /*ThinkingMsg*/

void ResizeWindow1 ( boolean reset )
{
  int w1, w2, h1, h2;

  SetWindow ( 1 );
  gwinxx = 110+(int)((float)(current_width-110)*fwinxx);
  gwinyy = (int)((float)(current_height)*fwinyy);

  w1 = gwinxx-1-110;  w2 = current_width-4-110-w1;
  h1 = gwinyy-1;      h2 = current_height-3-h1;
  edrect1[FIRST_GEOMW+0].x = 110;       edrect1[FIRST_GEOMW+0].y = 0;
  edrect1[FIRST_GEOMW+0].w = w1;        edrect1[FIRST_GEOMW+0].h = h1;
  edrect1[FIRST_GEOMW+1].x = 110+w1+3;  edrect1[FIRST_GEOMW+1].y = 0;
  edrect1[FIRST_GEOMW+1].w = w2;        edrect1[FIRST_GEOMW+1].h = h1;
  edrect1[FIRST_GEOMW+2].x = 110;       edrect1[FIRST_GEOMW+2].y = h1+3;
  edrect1[FIRST_GEOMW+2].w = w1;        edrect1[FIRST_GEOMW+2].h = h2;
  edrect1[FIRST_GEOMW+3].x = 110+w1+3;  edrect1[FIRST_GEOMW+3].y = h1+3;
  edrect1[FIRST_GEOMW+3].w = w2;        edrect1[FIRST_GEOMW+3].h = h2;
  edrect1[FIRST_GEOMW+4].x = 110;       edrect1[FIRST_GEOMW+4].y = 0;
  edrect1[FIRST_GEOMW+4].w = w1+w2+3;   edrect1[FIRST_GEOMW+4].h = h1+h2+3;
  edrect1[6].h = current_height;
  menurect1c[9].y = current_height-19;

  SetupParProj ( &edrect1[FIRST_GEOMW] );
  SetupPerspProj ( &edrect1[FIRST_GEOMW+3], reset );
  PictureIsOn = false;
  if ( RenderingIsOn )
    OffRendering ();
  ProjectSurfCP ();
} /*ResizeWindow1*/

void SwitchMenu1 ( int menu )
{
  boolean s;

  menu1 = menu;
  s = eswEdRendering;
  eswEdRendering = (menu1 == 2) &&
                   (eswSections || eswLines1 || eswLines2 || eswLight);
  if ( s != eswEdRendering )
    RedrawGeomWindows ();
  edrect1[6].redraw ( &edrect1[6] );
} /*SwitchMenu1*/

/* ///////////////////////////////////////////////////////////////////////// */
static Window the_window;

static void SendClientMessage ()
{
  XEvent event;

  event.type = ClientMessage;
  event.xclient.message_type = 0;
  event.xclient.format = 8;
  XSendEvent ( thedisplay, the_window, false, 0, &event );
} /*SendClientMessage*/

void OffRendering ()
{
  StopRendering ();
  SetWindow ( 1 );
  edrect1[6].redraw ( &edrect1[6] );
} /*OffRendering*/

void OnRendering ()
{
  struct tms start;

  SetWindow ( 1 );
  the_window = thewindow;
  PictureIsOn = false;
  edrect1[FIRST_GEOMW+3].redraw ( &edrect1[FIRST_GEOMW+3] );
  ResetRenderer ();
  BeginRendering ( &edrect1[FIRST_GEOMW+3] );
  edrect1[6].redraw ( &edrect1[6] );
  SendClientMessage ();
  ticks_per_second = sysconf ( _SC_CLK_TCK );
  lastupdate = times ( &start );
} /*OnRendering*/

void OnOffRendering ()
{
  if ( !RenderingIsOn )
    OnRendering ();
  else
    OffRendering ();
} /*OnOffRendering*/

void ProcessOtherEvent ()
{
  struct tms now;
  clock_t    clock1, clock2;

  if ( theevent.type == ClientMessage ) {
    if ( RenderingIsOn ) {
      PictureIsOn = true;
      clock1 = times ( &now );
      do {
        RenderLine ();
        clock2 = times ( &now );
      } while ( clock2-clock1 <= ticks_per_second/15 && RenderingIsOn );
      if ( clock2-lastupdate >= ticks_per_second || !RenderingIsOn ) {
        SetWindow ( 1 );
        edrect1[FIRST_GEOMW+3].redraw ( &edrect1[FIRST_GEOMW+3] );
        lastupdate = clock2;
      }
      if ( RenderingIsOn )
        SendClientMessage ();
      else
        edrect1[6].redraw ( &edrect1[6] );
    }
  }
} /*ProcessOtherEvent*/

