
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

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>

#include "pkgeom.h"
#include "xgedit.h"
#include "xgledit.h"

#include "gl2Dwin.h"


xge_widget *menu;
xge_2Dwinf wind;

point2f triangle[3] = {{-0.8,-0.8},{0.9,0.0},{-0.8,0.8}};
point2f saved[3];
boolean marking[3] = {false,false,false};
int current_point;

/* ///////////////////////////////////////////////////////////////////////// */
void RysujOkno ( xge_widget *er, boolean onscreen )
{
  int     i;

  glXWaitX ();

  xgle_DrawGeomWinBackground ( er, 0 );
  if ( wind.display_coord && wind.inside )
    xgle_2DwinfDrawCursorPos ( &wind, xge_xx, xge_yy );

  xgle_SetGLCameraf ( &wind.CPos );
  glBegin ( GL_POLYGON );
  glColor3f ( 1.0, 1.0, 1.0 );
  for ( i = 0; i < 3; i++ )
    glVertex2fv ( &triangle[i].x );
  glEnd ();

  glPointSize ( 3.0 );
  glBegin ( GL_POINTS );
  glColor3f ( 1.0, 1.0, 0.0 );
  for ( i = 0; i < 3; i++ )
    if ( !marking[i] )
      glVertex2fv ( &triangle[i].x );
  glColor3f ( 1.0, 0.3, 0.0 );
  for ( i = 0; i < 3; i++ )
    if ( marking[i] )
      glVertex2fv ( &triangle[i].x );
  glEnd ();

  glXWaitGL ();
  xgleCopyGLRect ( er->w, er->h, er->x, er->y );

  xge_DrawGeomWinSelectionRect ( er, &wind.selection_rect );
  xge_2DwinfDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*RysujOkno*/

/* ///////////////////////////////////////////////////////////////////////// */
int NearestPointFound ( CameraRecf *CPos, short x, short y )
{
  int     d, dmin;
  int     i;
  point2f q;

  dmin = xge_MINDIST+1;
  for ( i = 0; i < 3; i++ ) {
    CameraProjectPoint2f ( CPos, &triangle[i], &q );
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
  point2f q;

  SetPoint2f ( &q, (float)x, (float)y );
  CameraUnProjectPoint2f ( &wind.CPos, &q, &triangle[current_point] );
} /*MoveCurrentPoint*/

void SelectPoints ( CameraRecf *CPos, short x0, short x1, short y0, short y1 )
{
  int     i;
  point2f q;

  if ( x1-x0 >= 2 || y1-y0 >= 2 ) {
    for ( i = 0; i < 3; i++ ) {
      CameraProjectPoint2f ( CPos, &triangle[i], &q );
      if ( q.x >= x0 && q.x <= x1 && q.y >= y0 && q.y <= y1 )
        marking[i] = true;
    }
  }
  else {
    if ( NearestPointFound ( CPos, x0, y0 ) )
      marking[current_point] = true;
  }
} /*SelectPoints*/

void UnSelectPoints ( CameraRecf *CPos, short x0, short x1, short y0, short y1 )
{
  int     i;
  point2f q;

  if ( x1-x0 >= 2 || y1-y0 >= 2 ) {
    for ( i = 0; i < 3; i++ ) {
      CameraProjectPoint2f ( CPos, &triangle[i], &q );
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
      TransPoint2f ( &wind.gwtrans, &saved[i], &triangle[i] );
} /*TransformPoints*/

int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( er ) {
/*    printf ( "%d, 0x%x, %d, %d, %d\n", er->id, msg, key, x, y ); */
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
          xge_2DwinfEnableGeomWidget ( &wind, xge_2DWIN_PANNING_TOOL );
        else
          xge_2DwinfEnableGeomWidget ( &wind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        break;
    case 2:
        if ( wind.selecting_mode )
          xge_2DwinfEnableGeomWidget ( &wind, xge_2DWIN_SELECTING_TOOL );
        else
          xge_2DwinfEnableGeomWidget ( &wind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        break;
    case 4:
        if ( wind.moving_tool )
          xge_2DwinfEnableGeomWidget ( &wind, xge_2DWIN_MOVING_TOOL );
        else
          xge_2DwinfEnableGeomWidget ( &wind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        break;
    case 5:
        if ( wind.scaling_tool )
          xge_2DwinfEnableGeomWidget ( &wind, xge_2DWIN_SCALING_TOOL );
        else
          xge_2DwinfEnableGeomWidget ( &wind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        break;
    case 6:
        if ( wind.rotating_tool )
          xge_2DwinfEnableGeomWidget ( &wind, xge_2DWIN_ROTATING_TOOL );
        else
          xge_2DwinfEnableGeomWidget ( &wind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        break;
    case 7:
        if ( wind.shear_tool )
          xge_2DwinfEnableGeomWidget ( &wind, xge_2DWIN_SHEAR_TOOL );
        else
          xge_2DwinfEnableGeomWidget ( &wind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        break;
    default:
        break;
      }
      break;

  case xgemsg_2DWIN_PICK_POINT:
      return NearestPointFound ( &wind.CPos, x, y );

  case xgemsg_2DWIN_MOVE_POINT:
      MoveCurrentPoint ( er, x, y );
      break;

  case xgemsg_2DWIN_SELECT_POINTS:
      SelectPoints ( &wind.CPos,
                     wind.selection_rect.x0, wind.selection_rect.x1,
                     wind.selection_rect.y0, wind.selection_rect.y1 );
      break;

  case xgemsg_2DWIN_UNSELECT_POINTS:
      UnSelectPoints ( &wind.CPos,
                     wind.selection_rect.x0, wind.selection_rect.x1,
                     wind.selection_rect.y0, wind.selection_rect.y1 );
      break;

  case xgemsg_2DWIN_SAVE_POINTS:
      memcpy ( saved, triangle, 3*sizeof(point2f) );
      break;

  case xgemsg_2DWIN_TRANSFORM_POINTS:
      TransformPoints ();
      break;

  default:
      break;
    }
  }
  else {
/*    printf ( "0x%x, %d, %d, %d\n", msg, key, x, y ); */
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
      wind.er->msgproc ( wind.er, xgemsg_RESIZE, 0,
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
static char sw7text[] = "shear";

void init_edwin ( void )
{
  xge_widget *rp, *edr;
/*  int dw, dh, dwmm, dhmm; */

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
  rp = xge_NewSwitch ( 0, rp, 7, 109, 16, 0, 80, sw7text, &wind.shear_tool );
  menu = xge_NewMenu ( 0, NULL, 1, 110, xge_HEIGHT, 0, 0, rp );
  edr = xge_New2Dwinf ( 0, menu, 0, xge_WIDTH-110, xge_HEIGHT, 110, 0,
                        &wind, RysujOkno );
  xge_SetWinEdRect ( edr );
  CameraSetDepthRangef ( &wind.CPos, -1.0, 1.0 );
  xge_Redraw ();
} /*init_edwin*/

void destroy_edwin ( void )
{
} /*destroy_edwin*/

int main ( int argc, char *argv[] )
{
  xgle_Init ( argc, argv, CallBack, NULL,
              XGLE_WANTED, XGLE_NOT_NEEDED, XGLE_NOT_NEEDED );
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();  
  xgle_Cleanup ();  
  exit ( 0 );
} /*main*/   

