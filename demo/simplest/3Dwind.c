
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkgeom.h"
#include "xgedit.h"

#include "3Dwind.h"

xge_3Dwind wind;
xge_widget *menu;

point3d triangle[3] = {{-0.5,0.0,0.0},{0.5,-0.5,0.0},{0.7,0.7,0.2}};
point3d saved[3];
boolean marking[3] = {false, false, false};
int current_point;

/* ///////////////////////////////////////////////////////////////////////// */
void RysujOkno ( xge_widget *er, boolean onscreen )
{
  XPoint p[3];
  point3d q;
  int i;

  xge_DrawGeomWinBackground ( er );
  xge_3DwindDrawCursorPos ( &wind, er->id, xge_xx, xge_yy );

  for ( i = 0; i < 3; i++ ) {
    CameraProjectPoint3d ( &wind.CPos[er->id], &triangle[i], &q );
    p[i].x = (short)(q.x+0.5);
    p[i].y = (short)(q.y+0.5);
  }
  xgeSetForeground ( xgec_White );
  xgeFillPolygon ( Convex, 3, p );
  xgeSetForeground ( xgec_Yellow );
  for ( i = 0; i < 3; i++ )
    if ( !marking[i] )
      xgeFillRectangle ( 3, 3, p[i].x-1, p[i].y-1 );
  xgeSetForeground ( xgec_OrangeRed );
  for ( i = 0; i < 3; i++ )
    if ( marking[i] )
      xgeFillRectangle ( 3, 3, p[i].x-1, p[i].y-1 );

  xge_DrawGeomWinSelectionRect ( er, &wind.selection_rect );
  xge_3DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujOkno*/

/* ///////////////////////////////////////////////////////////////////////// */
int NearestPointFound ( CameraRecd *CPos, short x, short y )
{
  int     d, dmin;
  int     i;
  point3d q;

  dmin = xge_MINDIST+1;
  for ( i = 0; i < 3; i++ ) {
    CameraProjectPoint3d ( CPos, &triangle[i], &q );
    d = abs((int)(q.x+0.5)-x) + abs((int)(q.y+0.5)-y);
    if ( d < dmin ) {
      dmin = d;
      current_point = i;
    }
  }
  return dmin < xge_MINDIST;
} /*NearestPointFound*/

void MoveCurrentPoint ( xge_widget *er, short x, short y )
{
  point3d q;

  CameraProjectPoint3d ( &wind.CPos[er->id], &triangle[current_point], &q );
  q.x = (double)x;  q.y = (double)y;
  CameraUnProjectPoint3d ( &wind.CPos[er->id], &q, &triangle[current_point] );
} /*MoveCurrentPoint*/

void SelectPoints ( CameraRecd *CPos, short x0, short x1, short y0, short y1 )
{
  int     i;
  point3d q;

  if ( x1-x0 >= 2 || y1-y0 >= 2 ) {
    for ( i = 0; i < 3; i++ ) {
      CameraProjectPoint3d ( CPos, &triangle[i], &q );
      if ( q.x >= x0 && q.x <= x1 && q.y >= y0 && q.y <= y1 )
        marking[i] = true;
    }
  }
  else {
    if ( NearestPointFound ( CPos, x0, y0 ) )
      marking[current_point] = true;
  }
} /*SelectPoints*/

void UnSelectPoints ( CameraRecd *CPos, short x0, short x1, short y0, short y1 )
{
  int     i;
  point3d q;

  if ( x1-x0 >= 2 || y1-y0 >= 2 ) {
    for ( i = 0; i < 3; i++ ) {
      CameraProjectPoint3d ( CPos, &triangle[i], &q );
      if ( q.x >= x0 && q.x <= x1 && q.y >= y0 && q.y <= y1 )
        marking[i] = false;
    }
  }
  else {
    if ( NearestPointFound ( CPos, x0, y0 ) )
      marking[current_point] = false;
  }
} /*UnSelectPoints*/

void TransformPoints ( void )
{
  int i;

  for ( i = 0; i < 3; i++ )
    if ( marking[i] )
      TransPoint3d ( &wind.gwtrans, &saved[i], &triangle[i] );
} /*TransformPoints*/

int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( er ) {
/*    printf ( "%d, %d, %d, %d, %d\n", er->id, msg, key, x, y ); */
    switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
      switch ( er->id = 0 ) {
  case 0:
        xge_done = 1;
        break;

  default:
        break;
      }
      break;

case xgemsg_SWITCH_COMMAND:
      switch ( er->id ) {
  case 1:
        if ( wind.panning )
          xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_PANNING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case 2:
        if ( wind.selecting_mode )
          xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_SELECTING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case 4:
        if ( wind.moving_tool )
          xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_MOVING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case 5:
        if ( wind.scaling_tool )
          xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_SCALING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case 6:
        if ( wind.rotating_tool )
          xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_ROTATING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case 7:
        if ( wind.shear_tool )
          xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_SHEAR_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  default:
        break;
      }
      break;

case xgemsg_3DWIN_RESIZE:
case xgemsg_3DWIN_PROJCHANGE:
      /* nothing to do in this application, */
      /* another one might update the 3D object representation */
      break;

case xgemsg_3DWIN_PICK_POINT:
      return NearestPointFound ( &wind.CPos[er->id], x, y );

case xgemsg_3DWIN_MOVE_POINT:
      MoveCurrentPoint ( er, x, y );
      break;

case xgemsg_3DWIN_SELECT_POINTS:
      SelectPoints ( &wind.CPos[er->id],
                     wind.selection_rect.x0, wind.selection_rect.x1,
                     wind.selection_rect.y0, wind.selection_rect.y1 );
      break;

case xgemsg_3DWIN_UNSELECT_POINTS:
      UnSelectPoints ( &wind.CPos[er->id],
                       wind.selection_rect.x0, wind.selection_rect.x1,
                       wind.selection_rect.y0, wind.selection_rect.y1 );
      break;

case xgemsg_3DWIN_SAVE_POINTS:
      memcpy ( saved, triangle, 3*sizeof(point3d) );
      break;

case xgemsg_3DWIN_TRANSFORM_POINTS:
      TransformPoints ();
      break;

default:
      break;
    }
  }
  else {
/*    printf ( "%d, %d, %d, %d\n", msg, key, x, y ); */
    switch ( msg ) {
case xgemsg_KEY:
      switch ( key ) {
  case 'Q': case 'q':
        xge_done = 1;
        break;
  default:
        break;
      }
      break;

case xgemsg_RESIZE:
      wind.fww.er->msgproc ( wind.fww.er, xgemsg_RESIZE, 0,
                   (short)(xge_current_width-110), xge_current_height );
      menu->msgproc ( menu, xgemsg_RESIZE, 0, 110, xge_current_height );
      xge_Redraw ();
      break;

default:
      break;
    }
  }
  return 0;
} /*CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
static char b0[]      = "Quit";
static char sw1text[] = "pan/zoom";
static char sw2text[] = "select/unselect";
static char sw3text[] = "coordinates";
static char sw4text[] = "moving";
static char sw5text[] = "scaling";
static char sw6text[] = "rotating";

void init_edwin ( void )
{
  xge_widget *rp, *edr;
/*  int dw, dh, dwmm, dhmm; */

  pkv_InitScratchMem ( 65536 );
  rp = xge_NewButton ( 0, NULL, 0, 60, 18, 0, 0, b0 );
  rp = xge_NewSwitch ( 0, rp, 1, 109, 16, 0, xge_HEIGHT-16, sw1text, &wind.panning );
  xge_SetWidgetPositioning ( rp, 2, 0, -16 );
  rp = xge_NewSwitch ( 0, rp, 2, 109, 16, 0, xge_HEIGHT-36, sw2text, &wind.selecting_mode );
  xge_SetWidgetPositioning ( rp, 2, 0, -36 );
  rp = xge_NewSwitch ( 0, rp, 3, 109, 16, 0, xge_HEIGHT-56, sw3text, &wind.display_coord );
  xge_SetWidgetPositioning ( rp, 2, 0, -56 );
  rp = xge_NewSwitch ( 0, rp, 4, 109, 16, 0, 20, sw4text, &wind.moving_tool );
  rp = xge_NewSwitch ( 0, rp, 5, 109, 16, 0, 40, sw5text, &wind.scaling_tool );
  rp = xge_NewSwitch ( 0, rp, 6, 109, 16, 0, 60, sw6text, &wind.rotating_tool );
  menu = xge_NewMenu ( 0, NULL, 1, 110, xge_HEIGHT, 0, 0, rp );
  edr = xge_New3Dwind ( 0, menu, 0, xge_WIDTH-110, xge_HEIGHT, 110, 0,
                        &wind, RysujOkno, RysujOkno );
  xge_SetWinEdRect ( edr );
  xge_Redraw ();
/*
  dw = DisplayWidth ( xgedisplay, xgescreen );
  dh = DisplayHeight ( xgedisplay, xgescreen );
  dwmm = DisplayWidthMM ( xgedisplay, xgescreen );
  dhmm = DisplayHeightMM ( xgedisplay, xgescreen );
  printf ( "%d, %d, %d, %d\n", dw, dh, dwmm, dhmm );
*/
} /*init_edwin*/

void destroy_edwin ( void )
{
  pkv_DestroyScratchMem ();
} /*destroy_edwin*/

int main ( int argc, char *argv[] )
{
  xge_Init ( argc, argv, CallBack, NULL );
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();  
  xge_Cleanup ();  
  exit ( 0 );
} /*main*/   

