
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"

#include "xgedit.h"
#include "xgeprivate.h"


/* ///////////////////////////////////////////////////////////////////////// */
void xge_OrderSelectionRect ( Box2s *sel_rect )
{
  short a;

  if ( sel_rect->x0 > sel_rect->x1 ) {
    a = sel_rect->x0;
    sel_rect->x0 = sel_rect->x1;
    sel_rect->x1 = a;
  }
  if ( sel_rect->y0 > sel_rect->y1 ) {
    a = sel_rect->y0;
    sel_rect->y0 = sel_rect->y1;
    sel_rect->y1 = a;
  }
} /*xge_OrderSelectionRect*/

void xge_DrawGeomWinBackground ( xge_widget *er )
{
  unsigned long colour;

  switch ( er->state ) {
case xgestate_2DWIN_MOVINGPOINT:
case xgestate_2DWIN_ZOOMING:
case xgestate_2DWIN_PANNING:
case xgestate_2DWIN_SELECTING:
case xgestate_2DWIN_UNSELECTING:
case xgestate_2DWIN_TGSELECTING:
case xgestate_2DWIN_MOVING_GEOM_WIDGET:
case xgestate_2DWIN_USING_GEOM_WIDGET:
case xgestate_2DWIN_ALTUSING_GEOM_WIDGET:
case xgestate_3DWIN_MOVINGPOINT:
case xgestate_3DWIN_ZOOMING:
case xgestate_3DWIN_PARPANNING:
case xgestate_3DWIN_PARZOOMING:
case xgestate_3DWIN_TURNING_VIEWER:
case xgestate_3DWIN_SELECTING:
case xgestate_3DWIN_UNSELECTING:
case xgestate_3DWIN_TGSELECTING:
case xgestate_3DWIN_MOVING_GEOM_WIDGET:
case xgestate_3DWIN_USING_GEOM_WIDGET:
case xgestate_3DWIN_ALTUSING_GEOM_WIDGET:
case xgestate_KNOTWIN_MOVINGKNOT:
case xgestate_KNOTWIN_PANNING:
case xgestate_KNOTWIN_ZOOMING:
case xgestate_T2KNOTWIN_MOVINGKNOT_U:
case xgestate_T2KNOTWIN_MOVINGKNOT_V:
case xgestate_T2KNOTWIN_MOVING_POINT:
case xgestate_T2KNOTWIN_PANNING:
case xgestate_T2KNOTWIN_ZOOMING:
case xgestate_T2KNOTWIN_SELECTING:
case xgestate_T2KNOTWIN_UNSELECTING:
case xgestate_T2KNOTWIN_TGSELECTING:
    colour = xgec_GEOM_WIN_BACK_A;
    break;
default:
    colour = xgec_GEOM_WIN_BACK_B;
    break;
  }
  xgeSetForeground ( colour );
  xgeFillRectangle ( er->w-2, er->h-2, er->x+1, er->y+1 );
} /*xge_DrawGeomWinBackground*/

void xge_DrawGeomWinFrame ( xge_widget *er, boolean onscreen )
{
  xgeSetForeground ( xgec_GEOM_WIN_FRAME );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawGeomWinFrame*/

void xge_DrawGeomWinSelectionRect ( xge_widget *er, Box2s *sel_rect )
{
  short         x0, x1, y0, y1;
  static char   dla[] = {2,2};
  static char   dlb[] = {2,4};
  unsigned long colour;

  if ( sel_rect->x0 > sel_rect->x1 )
    { x0 = sel_rect->x1;  x1 = sel_rect->x0; }
  else
    { x0 = sel_rect->x0;  x1 = sel_rect->x1; }
  if ( sel_rect->y0 > sel_rect->y1 )
    { y0 = sel_rect->y1;  y1 = sel_rect->y0; }
  else
    { y0 = sel_rect->y0;  y1 = sel_rect->y1; }
  if ( x1-x0 < 1 && y1-y0 < 1 )  /* do not draw too small rectangles */
    return;
  switch ( er->state ) {
case xgestate_2DWIN_SELECTING:
case xgestate_3DWIN_SELECTING:
case xgestate_T2KNOTWIN_SELECTING:
    colour = xgec_GEOM_WIN_SEL_A;
    goto draw_it;
case xgestate_2DWIN_UNSELECTING:
case xgestate_3DWIN_UNSELECTING:
case xgestate_T2KNOTWIN_UNSELECTING:
    colour = xgec_GEOM_WIN_SEL_B;
draw_it:
    xgeSetForeground ( colour );
    xgeSetDashes ( 2, dla, 0 );
    xgeSetLineAttributes ( 1, LineOnOffDash, CapRound, JoinRound );
    xgeDrawRectangle ( x1-x0+1, y1-y0+1, x0, y0 );
    xgeSetLineAttributes ( 1, LineSolid, CapRound, JoinRound );
    break;
case xgestate_2DWIN_TGSELECTING:
case xgestate_3DWIN_TGSELECTING:
case xgestate_T2KNOTWIN_TGSELECTING:
    xgeSetLineAttributes ( 1, LineOnOffDash, CapRound, JoinRound );
    xgeSetForeground ( xgec_GEOM_WIN_SEL_A );
    xgeSetDashes ( 2, dlb, 0 );
    xgeDrawRectangle ( x1-x0+1, y1-y0+1, x0, y0 );
    xgeSetForeground ( xgec_GEOM_WIN_SEL_B );
    xgeSetDashes ( 2, dlb, 2 );
    xgeDrawRectangle ( x1-x0+1, y1-y0+1, x0, y0 );
    xgeSetLineAttributes ( 1, LineSolid, CapRound, JoinRound );
    break;
default:
    break;
  }
} /*xge_DrawGeomWinSelectionRect*/

