
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
#include <X11/cursorfont.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "xgedit.h"
#include "bookg1holef.h"

#include "datagen.h"
#include "drawg1hole.h"
#include "edg1hole.h"

xge_widget *menu0, *menu1;
xge_widget *menu00list, *menu01list, *menu02list;

xge_int_widget wdg_holek;

/* ///////////////////////////////////////////////////////////////////////// */
void RysujOkno ( xge_widget *er, boolean onscreen )
{
  int id;

  id = er->id;
  xge_DrawGeomWinBackground ( er );

  if ( id < 3 )
    xge_3DwinfDrawCursorPos ( &swind, id, xge_xx, xge_yy );

  DrawSurface ( id );

  xge_DrawGeomWinSelectionRect ( er, &swind.selection_rect );
  if ( id < 3 )
    xge_3DwinfDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujOkno*/

/* ///////////////////////////////////////////////////////////////////////// */
static char *InfoMsg[7] =
  { "Filling polygonal holes with tangent plane continuity",
    "using the procedure described in the book (in Polish)",
    "Podstawy modelowania krzywych i powierzchni"
    "",
    "This program is a part of the BSTools package.",
    "(C) Copyright by Przemyslaw Kiciak, 2005,2007",
    NULL };

static char btn0[] = "Quit";
static char btn1[] = "About";
static char btn2[] = "Data";
static char btn3[] = "Edit";
static char btn4[] = "View";
static char mtext2[] = "beta";

static char intw0text[] = "number of sides";
static char sw0text[]   = "select/unselect";
static char sw1text[]   = "move";
static char sw2text[]   = "scale";
static char sw3text[]   = "rotate";
static char sw20text[]  = "shear";
static char swtext11[]  = "auto centre";
static char swtext12[]  = "B-spline net";
static char swtext13[]  = "patches";
static char swtext14[]  = "aux. curves";
static char swtext15[]  = "final curves";
static char swtext16[]  = "aux. patches";
static char swtext17[]  = "final patches";
static char sw18text[]  = "pan & zoom";
static char sw19text[]  = "coordinates";

void init_edwin ( void )
{
  xge_widget *w1, *w2;

  pkv_InitScratchMem ( SCRATCH_MEM_SIZE );

    /* setup the four geometry windows */
  w1 = xge_New3Dwinf ( 0, NULL, 0, xge_WIDTH-110, xge_HEIGHT-20, 110, 20,
                       &swind, RysujOkno, RysujOkno );
  xge_3DwinfSetDefBBox ( &swind, -2.0, 2.0, -2.0, 2.0, -2.0, 2.0 );
  xge_3DwinfInitProjections ( &swind, -2.0, 2.0, -2.0, 2.0, -2.0, 2.0 );

    /* setup the top menu */
  w2 = xge_NewButton ( 0, NULL, 0, 60, 19, 0, 0, btn0 );
  w2 = xge_NewButton ( 0, w2, 1, 60, 19, (short)(xge_WIDTH-60), 0, btn1 );
  xge_SetWidgetPositioning ( w2, 1, -60, 0 );
  w2 = xge_NewButton ( 0, w2, 2, 60, 19,  62, 0, btn2 );
  w2 = xge_NewButton ( 0, w2, 3, 60, 19, 124, 0, btn3 );
  w2 = xge_NewButton ( 0, w2, 4, 60, 19, 186, 0, btn4 );
  menu1 = xge_NewMenu ( 0, w1, 7, xge_WIDTH, 20, 0, 0, w2 );

    /* setup the side menu 00 - Data */
  w2 = xge_NewIntWidget ( 0, NULL, 0, 109, 19, 0, 20, 2, MAX_K+1,
                          &wdg_holek, intw0text, &hole_k );
  w2 = xge_NewSlidebarf ( 0, w2,  6, 109, 10, 0,  48, &dataparam[0] );
  w2 = xge_NewSlidebarf ( 0, w2,  7, 109, 10, 0,  66, &dataparam[1] );
  w2 = xge_NewSlidebarf ( 0, w2,  8, 109, 10, 0,  84, &dataparam[2] );
  w2 = xge_NewSwitch ( 0, w2, 18, 109, 16, 0, xge_HEIGHT-36,
                      sw18text, &swind.panning );
  xge_SetWidgetPositioning ( w2, 2, 0, -36 );
  w2 = xge_NewSwitch ( 0, w2, 19, 109, 16, 0, xge_HEIGHT-16,
                      sw19text, &swind.display_coord );
  xge_SetWidgetPositioning ( w2, 2, 0, -16 );
  menu00list = w2;

    /* setup the side menu 01 - Edit */
  w2 = xge_NewSwitch ( 0, NULL, 0, 109, 16, 0, 20, sw0text, &swind.selecting_mode );
  w2 = xge_NewSwitch ( 0, w2, 1, 109, 16, 0, 40, sw1text, &swind.moving_tool );
  w2 = xge_NewSwitch ( 0, w2, 2, 109, 16, 0, 60, sw2text, &swind.scaling_tool );
  w2 = xge_NewSwitch ( 0, w2, 3, 109, 16, 0, 80, sw3text, &swind.rotating_tool );
  w2 = xge_NewSwitch ( 0, w2, 20, 109, 16, 0, 100, sw20text, &swind.shear_tool );
  w2 = xge_NewTextWidget ( 0, w2, 0, 109, 16, 0, 116, mtext2 );
  w2 = xge_NewSlidebarf ( 0, w2,  9, 109, 10, 0, 132, &beta1 );
  w2 = xge_NewSlidebarf ( 0, w2, 10, 109, 10, 0, 150, &beta2 );
  w2 = xge_NewSwitch ( 0, w2, 11, 109, 16, 0, 168, swtext11, &autocpoint );
  w2 = xge_NewSwitch ( 0, w2, 18, 109, 16, 0, xge_HEIGHT-36,
                      sw18text, &swind.panning );
  xge_SetWidgetPositioning ( w2, 2, 0, -36 );
  w2 = xge_NewSwitch ( 0, w2, 19, 109, 16, 0, xge_HEIGHT-16,
                      sw19text, &swind.display_coord );
  xge_SetWidgetPositioning ( w2, 2, 0, -16 );
  menu01list = w2;
    /* setup the side menu 02 - View */
  w2 = xge_NewSwitch ( 0, NULL, 12, 109, 16, 0, 20, swtext12, &bsnet );
  w2 = xge_NewSwitch ( 0, w2, 13, 109, 16, 0, 40, swtext13, &bezp );
  w2 = xge_NewSwitch ( 0, w2, 14, 109, 16, 0, 60, swtext14, &auxcurv );
  w2 = xge_NewSwitch ( 0, w2, 15, 109, 16, 0, 80, swtext15, &starcurv );
  w2 = xge_NewSwitch ( 0, w2, 16, 109, 16, 0, 100, swtext16, &auxpatch );
  w2 = xge_NewSwitch ( 0, w2, 17, 109, 16, 0, 120, swtext17, &targpatch );
  w2 = xge_NewSwitch ( 0, w2, 18, 109, 16, 0, xge_HEIGHT-36,
                      sw18text, &swind.panning );
  xge_SetWidgetPositioning ( w2, 2, 0, -36 );
  w2 = xge_NewSwitch ( 0, w2, 19, 109, 16, 0, xge_HEIGHT-16,
                      sw19text, &swind.display_coord );
  xge_SetWidgetPositioning ( w2, 2, 0, -16 );
  menu02list = w2;
  menu0 = xge_NewMenu ( 0, menu1, 6, 110, xge_HEIGHT-20, 0, 20, menu00list );
  
  ResetObject ();
  xge_SetWinEdRect ( menu0 );
  xge_Redraw ();
  xge_DisplayInfoMessage ( InfoMsg, -1 );
} /*init_edwin*/

void destroy_edwin ( void )
{
  printf ( "Scratch memory used: %d out of %d bytes\n",
           (int)pkv_MaxScratchTaken(), SCRATCH_MEM_SIZE );
  pkv_DestroyScratchMem ();
} /*destroy_edwin*/

/* ///////////////////////////////////////////////////////////////////////// */
int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( er ) {
    switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
      switch ( er->id ) {
  case 0:  /* Quit */
        xge_done = true;
        break;
  case 1:  /* About */
        xge_DisplayInfoMessage ( InfoMsg, -1 );
        break;
  case 2:  /* Data */
        xge_SetMenuWidgets ( menu0, menu00list, true );
        break;
  case 3:  /* Edit */
        xge_SetMenuWidgets ( menu0, menu01list, true );
        break;
  case 4:  /* View */
        xge_SetMenuWidgets ( menu0, menu02list, true );
        break;
 default:
        return 0;
      }
      break;

case xgemsg_SLIDEBAR_COMMAND:
      switch ( er->id ) {
  case 6:  case 7:  case 8:  /* data generator parameters */
        InitHole ( hole_k, dataparam[0], dataparam[1], dataparam[2] );
        UpdateSurface ();
        ResizeObject ();
        xge_SetClipping ( swind.fww.er );
        swind.fww.er->redraw ( swind.fww.er, true );
        break;

  case 9:  case 10:  /* construction parameters */
        UpdateSurface ();
        xge_SetClipping ( swind.fww.er );
        swind.fww.er->redraw ( swind.fww.er, true );
        break;

  default:
        break;
      }
      break;

case xgemsg_INT_WIDGET_COMMAND:
      if ( er->id == 0 ) {
        switch ( key ) {
    case 2:
          xge_DisplayErrorMessage ( "Less than 3 is not allowed.", -1 );
          return 0;
    case 3:
setup3:
          UstawK3 ();
          break;
    case 4:
          if ( key < hole_k ) goto setup3;
    case 5:
          UstawK5 ();
          break;
    case 6:
setup6:
          UstawK6 ();
          break;
    case 7:
          if ( key < hole_k ) goto setup6;
    case 8:
          UstawK8 ();
          break;
    case MAX_K+1:  /* 9 */
          xge_DisplayErrorMessage ( "More than 8 is not allowed.", -1 );
          return 0;
    default:
          return 0;
        }
        xge_Redraw ();
      }
      break;

case xgemsg_SWITCH_COMMAND:
      switch ( er->id ) {
  case 0:  /* select/unselect */
        if ( swind.selecting_mode )
          xge_3DwinfEnableGeomWidget ( &swind, xge_3DWIN_SELECTING_TOOL );
        else
          xge_3DwinfEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case 1:  /* move */
        if ( swind.moving_tool )
          xge_3DwinfEnableGeomWidget ( &swind, xge_3DWIN_MOVING_TOOL );
        else
          xge_3DwinfEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case 2:  /* scale */
        if ( swind.scaling_tool )
          xge_3DwinfEnableGeomWidget ( &swind, xge_3DWIN_SCALING_TOOL );
        else
          xge_3DwinfEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case 3:  /* rotate */
        if ( swind.rotating_tool )
          xge_3DwinfEnableGeomWidget ( &swind, xge_3DWIN_ROTATING_TOOL );
        else
          xge_3DwinfEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case 20:  /* shear */
        if ( swind.shear_tool )
          xge_3DwinfEnableGeomWidget ( &swind, xge_3DWIN_SHEAR_TOOL );
        else
          xge_3DwinfEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case 11:
        UstawAutoCPoint ();
        swind.fww.er->redraw ( swind.fww.er, true );
        break;
  case 12:  case 13:  case 14:  case 15:  case 16:  case 17:
        xge_Redraw ();
        break;
  case 18:  /*panning */
        if ( swind.panning )
          xge_3DwinfEnableGeomWidget ( &swind, xge_3DWIN_PANNING_TOOL );
        else
          xge_3DwinfEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case 19:  /* coordinates */
        break;
  default:
        return 0;
      }
      break;

case xgemsg_3DWIN_RESIZE:
case xgemsg_3DWIN_PROJCHANGE:
      ResizeObject ();
      break;

case xgemsg_3DWIN_PICK_POINT:
      return FindNearestPoint ( er->id, x, y, xge_MINDIST );

case xgemsg_3DWIN_MOVE_POINT:
      SetCPoint ( er->id, x, y );
      break;

case xgemsg_3DWIN_SELECT_POINTS:
      SelectCPoints ( er->id,
                      swind.selection_rect.x0, swind.selection_rect.x1,
                      swind.selection_rect.y0, swind.selection_rect.y1 );
      break;

case xgemsg_3DWIN_UNSELECT_POINTS:
      UnSelectCPoints ( er->id,
                        swind.selection_rect.x0, swind.selection_rect.x1,
                        swind.selection_rect.y0, swind.selection_rect.y1 );
      break;

case xgemsg_3DWIN_SAVE_POINTS:
      SaveCPoints ();
      break;

case xgemsg_3DWIN_TRANSFORM_POINTS:
      TransformCPoints ( &swind.gwtrans );
      break;

case xgemsg_3DWIN_FIND_REFBBOX:
      FindBoundingBox ( &swind.RefBBox );
      break;

default:
      break;
    }
  }
  else {
    switch ( msg ) {
case xgemsg_RESIZE:
      swind.fww.er->msgproc ( swind.fww.er, xgemsg_RESIZE, 0,
                              (short)(xge_current_width-110),
                              (short)(xge_current_height-20) );
      menu0->msgproc ( menu0, xgemsg_RESIZE, 0,
                110, (short)(xge_current_height-20) );
      menu1->msgproc ( menu1, xgemsg_RESIZE, 0, xge_current_width, 20 );
      ResizeObject ();
      xge_Redraw ();
      break;

case xgemsg_KEY:
      switch ( key ) {
  case 'M':
        xgeResizeWindow ( xge_MAX_WIDTH, xge_MAX_HEIGHT );
        break;
  case 'm':
        xgeResizeWindow ( xge_WIDTH, xge_HEIGHT );
        break;
  case 'D': case 'd':
        xgeMoveWindow ( 0, 512 );
        break;
  case 'Q': case 'q':
        xge_done = 1;
        break;
  default:
        break;
      }
      break;

default:
      break;
    }
  }
  return 0;
} /*CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
int main ( int argc, char *argv[] )
{
  xge_Init ( argc, argv, CallBack, NULL );
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();
  xge_Cleanup ();
  exit ( 0 );
} /*main*/

