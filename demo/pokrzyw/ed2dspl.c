
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "xgedit.h"
#include "bsfile.h"

#include "spline2d.h"
#include "ed2dwidgets.h"
#include "ed2dspl.h"

xge_widget *menu0, *menu1, *popup00, *popup01, *popup02, *popup03;
xge_widget *menu00list, *menu01list;

xge_int_widget ddegree, ccdens;

xge_listbox filelist, dirlist;
xge_string_ed filename_editor;
const char file_filter[] = "*.bs";
const char file_ext[] = ".bs";
char current_directory[MAX_PATH_LGT+1] = "";
char current_dir[MAX_PATH_SHRT+1] = "";
char filename[MAX_FILENAME_LGT+1] = "";


/* additional states */
#define STATE_MOVINGKNOT   (xgestate_LAST+1)
#define STATE_ZOOMINGKNOTS (xgestate_LAST+2)
#define STATE_PANNINGKNOTS (xgestate_LAST+3)

#define MIN_PARZOOM  0.1
#define MAX_PARZOOM 10.0

/* ///////////////////////////////////////////////////////////////////////// */
void RysujOkno ( xge_widget *er, boolean onscreen )
{
  xge_2Dwind *_2Dwin;

  _2Dwin = er->data0;
  xge_DrawGeomWinBackground ( er );
  if ( _2Dwin->inside && _2Dwin->display_coord )
    xge_2DwindDrawCursorPos ( _2Dwin, xge_xx, xge_yy );
  if ( convh )
    DisplayConvh ();
  if ( baza )
    DisplayBasis ( er->y+er->h-11, 20-er->h );
  if ( lamana )
    DisplayCPolygon ();
  if ( bezpoly )
    DisplayBezPoly ();
  if ( lamana )
    DisplayCPoints ();
  if ( polarf )
    DisplayPolarf ();
  if ( nurbs )
    DisplayAuxPoints ();
  if ( curvgr )
    DisplayCurvGraph ();
  if ( curve )
    DisplayCurve ();
  xge_DrawGeomWinSelectionRect ( er, &_2Dwin->selection_rect );
  xge_2DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujOkno*/

void ResetCurveMapping ( void )
{
  xge_2DwindInitProjection ( &cwind, -1.0, 1.0, -1.0, 1.0 );
  UstawFunkcje ();
  ProjectCurve ();
} /*ResetCurveMapping*/

void FindCurveMapping ( void )
{
  Box2d box;

  FindBoundingBox ( &box );
  xge_2DwindInitProjection ( &cwind, box.x0, box.x1, box.y0, box.y1 );
  UstawFunkcje ();
  ProjectCurve ();
} /*FindCurveMapping*/

void RedrawCurveAndKnotsArea ()
{
  xge_SetClipping ( cwind.er );
  cwind.er->redraw ( cwind.er, true );
  xge_SetClipping ( kwind.er );
  kwind.er->redraw ( kwind.er, true );
} /*RedrawCurveAndKnotsArea*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean FilenameCorrect ( const char *filename )
{
  return (boolean)(filename[0] != 0);
} /*FilenameCorrect*/

boolean SaveBSCurve ( char *filename )
{
  void   *sp;
  int    lfn, lex;
  double *buf;

  sp = pkv_GetScratchMemTop ();
  lfn = strlen ( filename );
  lex = strlen ( file_ext );
  if ( strcmp ( file_ext, &filename[lfn-lex] ) )
    strcpy ( &filename[lfn], file_ext );
  if ( bsf_OpenOutputFile ( filename, false ) ) {
    if ( nurbs ) {
try_nurbs:
      bsf_WriteBSplineCurved ( 2, 3, true, kwind.degree, kwind.lastknot,
                               knots, kwind.closed, &cpoints[0].x,
                               cpmark, NULL );
    }
    else {
      buf = pkv_GetScratchMem ( (kwind.lastknot-kwind.degree)*sizeof(point2d) );
      if ( !buf )
        goto try_nurbs;
      pkv_Selectd ( kwind.lastknot-kwind.degree, 2, 3, 2,
                    &cpoints[0].x, buf );
      bsf_WriteBSplineCurved ( 2, 2, false, kwind.degree, kwind.lastknot,
                               knots, kwind.closed, buf, cpmark, NULL );
    }
    bsf_CloseOutputFile ();
    pkv_SetScratchMemTop ( sp );
    return true;
  }
  else {
    pkv_SetScratchMemTop ( sp );
    return false;
  }
} /*SaveBSCurve*/

boolean ReadBSCurve ( char *filename )
{
  void    *sp;
  point4d *cp;
  double  *kn;
  int     maxcp, spdimen, degree, lastknot, i;
  boolean result, closed, rational;
  byte    *mkcp;

  sp = pkv_GetScratchMemTop ();
  maxcp = MAX_DEGREE*MAX_KNOTS;
  cp = pkv_GetScratchMem ( maxcp*sizeof(point4d) );
  kn = pkv_GetScratchMemd ( MAX_KNOTS );
  mkcp = pkv_GetScratchMem ( maxcp );
  if ( !cp || !kn )
    goto failure;

  if ( bsf_OpenInputFile ( filename ) ) {
    result = bsf_ReadBSplineCurve4d ( MAX_DEGREE, MAX_KNOTS, maxcp,
               &degree, &lastknot, kn, &closed, cp, &spdimen, &rational,
               mkcp, NULL );
    bsf_CloseInputFile ();
    if ( result ) {
      kwind.lastknot = lastknot;
      memcpy ( knots, kn, (lastknot+1)*sizeof(double) );
      kwind.degree = degree;
      npoints = lastknot-degree;
      memcpy ( cpmark, mkcp, npoints );
      pkv_Selectd ( npoints, 3, 4, 3, cp, cpoints );
      if ( spdimen <= 2 ) {
        for ( i = 0; i < npoints; i++ )
          cpoints[i].z = 1.0;
      }
      else if ( spdimen == 3 || spdimen == 4 ) {
        for ( i = 0; i < npoints; i++ ) {
          cpoints[i].z = cp[i].w;
          if ( !cpoints[i].z )
            goto failure;
        }
      }
      kwind.closed = closed;
      nurbs = rational;
      if ( closed ) {
        kwind.clcK = lastknot - 2*degree;
        kwind.clcT = knots[lastknot-degree] - knots[degree];
      } 
    }
    else
      bsf_PrintErrorLocation ();
  }
  else
    goto failure;

  funkcja = false;
  pkv_SetScratchMemTop ( sp );
  return result;

failure:
  ResetObject ();
  pkv_SetScratchMemTop ( sp );
  return false;
} /*ReadBSCurve*/

void SetCurrentDir ( void )
{
  int l;

  getcwd ( current_directory, MAX_PATH_LGT+1 );
  l = strlen ( current_directory );
  if ( l <= MAX_PATH_SHRT )
    strcpy ( current_dir, current_directory );
  else {
    strcpy ( current_dir, "..." );
    strcpy ( &current_dir[3], &current_directory[l-MAX_PATH_SHRT+3] );
  }
} /*SetCurrentDir*/

/* ///////////////////////////////////////////////////////////////////////// */
void init_edwin ( void )
{
  xge_widget *w;

  pkv_InitScratchMem ( SCRATCH_MEM_SIZE );

  /* setup the menu 00 widgets */
  w = xge_NewIntWidget ( 0, NULL, intw00DEG, 109, 19, 0, 20, 0, MAX_DEGREE+1,
                         &ddegree, intw2text, &kwind.degree );
  w = xge_NewSwitch ( 0, w, sw00SELECT, 109, 16, 0,  40,
                      sw0text, &cwind.selecting_mode );
  w = xge_NewSwitch ( 0, w, sw00MOVE, 109, 16, 0,  60,
                      sw1text, &cwind.moving_tool );
  w = xge_NewSwitch ( 0, w, sw00SCALE, 109, 16, 0,  80,
                      sw2text, &cwind.scaling_tool );
  w = xge_NewSwitch ( 0, w, sw00ROTATE, 109, 16, 0, 100,
                      sw3text, &cwind.rotating_tool );
  w = xge_NewSwitch ( 0, w, sw00SHEAR, 109, 16, 0, 120,
                      txtShear, &cwind.shear_tool );
  w = xge_NewSwitch ( 0, w, sw00NURBS, 109, 16, 0, 140,
                      sw9text, &nurbs );
  w = xge_NewSwitch ( 0, w, sw00CLOSED, 109, 16, 0, 160,
                      sw10text, &kwind.closed );
  w = xge_NewSwitch ( 0, w, sw00FUNCTION, 109, 16, 0, 180,
                      sw11text, &funkcja );
  w = xge_NewSwitch ( 0, w, sw00UNIFORM, 109, 16, 0, 200,
                      sw19text, &uniform );
  w = xge_NewButton ( 0, w, btn00REFINE, 60, 19, 0, 220, b4text );
  w = xge_NewButton ( 0, w, btn00RESET, 60, 19, 0, 240, b5text );
    /* these 3 switches are to be moved after window resizing */
  w = xge_NewSwitch ( 0, w, sw00PANZOOM, 109, 16, 0, xge_HEIGHT-56,
                           sw16text, &cwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( 0, w, sw00COORD, 109, 16, 0, xge_HEIGHT-36,
                           sw17text, &cwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( 0, w, sw00MOVEMANY , 109, 16, 0, xge_HEIGHT-16,
                           sw18text, &kwind.moving_many );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );

  menu00list = w;

  /* setup the menu 01 widgets */
  w = xge_NewSwitch ( 0, NULL, sw01CURVE, 109, 16, 0,  20, sw4text, &curve );
  w = xge_NewSwitch ( 0, w, sw01POLYLINE, 109, 16, 0,  40, sw5text, &lamana );
  w = xge_NewSwitch ( 0, w, sw01BEZPOLY, 109, 16, 0,  60, sw6text, &bezpoly );
  w = xge_NewSwitch ( 0, w, sw01TICKS, 109, 16, 0,  80, sw7text, &ticks );
  w = xge_NewSwitch ( 0, w, sw01CONVH, 109, 16, 0, 100, sw8text, &convh );
  w = xge_NewSwitch ( 0, w, sw01BASIS, 109, 16, 0, 120, sw12text, &baza );
  w = xge_NewSwitch ( 0, w, sw01POLARF, 109, 16, 0, 140, sw13text, &polarf );
  w = xge_NewSwitch ( 0, w, sw01CURVGR, 109, 16, 0, 160, sw14text, &curvgr );
  w = xge_NewSlidebard ( 0, w, sl01CURVGRSC, 109, 10, 0, 180, &curvgraphscale );
  w = xge_NewIntWidget ( 0, w, intw01CURVGRDENS, 109, 19, 0, 194, 1, 32,
                         &ccdens, intw15text, &curv_graph_dens );
    /* these 3 switches are to be moved after window resizing */
  w = xge_NewSwitch ( 0, w, sw00PANZOOM, 109, 16, 0, xge_HEIGHT-56,
                           sw16text, &cwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( 0, w, sw00COORD, 109, 16, 0, xge_HEIGHT-36,
                           sw17text, &cwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( 0, w, sw00MOVEMANY, 109, 16, 0, xge_HEIGHT-16,
                           sw18text, &kwind.moving_many );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  /* setup the windows */
    /* the menu with the above widgets */
  menu01list = w;
  menu0 = xge_NewMenu ( 0, NULL, MENU00, 110, xge_HEIGHT-20, 0, 20, menu00list );

  w = xge_NewButton ( 0, NULL, btn02FILE, 60, 19,  0,  0, txtFile );
  w = xge_NewButton ( 0, w, btn02EDIT, 60, 19, 62,  0, b1text );
  w = xge_NewButton ( 0, w, btn02VIEW, 60, 19, 124, 0, b2text );
  w = xge_NewButton ( 0, w, btn02ABOUT, 60, 19, xge_WIDTH-60, 0, b3text );
  xge_SetWidgetPositioning ( w, 1, -60, 0 );
    /* the menu with the above widgets */
  menu1 = xge_NewMenu ( 0, menu0, MENU02, xge_WIDTH, 20, 0, 0, w );

  w = xge_NewButton ( 0, NULL, btn00pOPEN, 76, 19, 2, 22, txtOpen );
  w = xge_NewButton ( 0, w, btn00pSAVE, 76, 19, 2, 42, txtSave );
  w = xge_NewButton ( 0, w, btn00pSAVEAS, 76, 19, 2, 62, txtSaveAs );
  w = xge_NewButton ( 0, w, btn00pEXPORT, 76, 19, 2, 82, txtExport );
  w = xge_NewButton ( 0, w, btn00pEXIT, 76, 19, 2, 102, txtExit );
  popup00 = xge_NewFMenu ( 0, NULL, POPUP00, 80, 103, 0, 20, w );
  popup00->msgproc = xge_PopupMenuMsg;

  /* setup the Open file menu */
  w = xge_NewTextWidget ( 0, NULL, 0, 380, 16, 20+10, 40+10, current_dir );
  w = xge_NewButton ( 0, w, btn01pOPEN, 76, 19, 20+82, 40+180-30, txtOpen );
  w = xge_NewButton ( 0, w, btn01pCANCEL, 76, 19, 20+242, 40+180-30, txtCancel );
  w = xge_NewListBox ( 0, w, lb01pDIRLIST, 180, 99, 20+10, 40+40, &dirlist );
  w = xge_NewListBox ( 0, w, lb01pFILELIST, 180, 99, 220+10, 40+40, &filelist );
  popup01 = xge_NewFMenu ( 0, NULL, POPUP01, 400, 180, 20, 40, w );

  /* setup the Save as file menu */
  w = xge_NewTextWidget ( 0, NULL, 0, 380, 16, 20+10, 40+10, current_dir );
  w = xge_NewButton ( 0, w, btn02pSAVE, 76, 19, 20+82, 40+180-30, txtSave );
  w = xge_NewButton ( 0, w, btn02pCANCEL, 76, 19, 20+242, 40+180-30, txtCancel );
  w = xge_NewListBox ( 0, w, lb02pDIRLIST, 180, 99, 20+10, 40+40, &dirlist );
  w = xge_NewTextWidget ( 0, w, txt02pSAVE_AS, 120, 16, 220+10, 40+38, txtSaveAs );
  w = xge_NewStringEd ( 0, w, txted02pFILENAME, 180, 19, 220+10, 40+56,
                        MAX_FILENAME_LGT, filename, &filename_editor );
  popup02 = xge_NewFMenu ( 0, NULL, POPUP02, 400, 180, 20, 40, w );

  /* setup the Exit menu */
  w = xge_NewTextWidget ( 0, NULL, 0, 240, 16, 110, 70, MsgReconsider  );
  w = xge_NewButton ( 0, w, btn03EXIT, 58, 19, 110, 100, txtExit );
  w = xge_NewButton ( 0, w, btn03SAVE, 58, 19, 210, 100, txtSave );
  w = xge_NewButton ( 0, w, btn03CANCEL, 58, 19, 310, 100, txtCancel );
  popup03 = xge_NewFMenu ( 0, NULL, POPUP03, 300, 70, 90, 60, w );
  popup03->msgproc = xge_PopupMenuMsg;

    /* the knot window */
  kwind.er = xge_NewKnotWind ( 0, menu1, KNOTWIN, xge_WIDTH-110, 56,
                               110, xge_HEIGHT-56, &kwind, MAX_KNOTS, knots );
    /* the window with the curve */
  cwind.er = xge_New2Dwind ( 0, kwind.er, CURVEWIN, xge_WIDTH-110,
                             xge_HEIGHT-78, 110, 20, &cwind, RysujOkno );
  xge_SetWinEdRect ( cwind.er );

  ResetObject ();
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
  boolean uu;

  if ( er ) {
    switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
      switch ( er->id ) {
  case btn02FILE:  /* File */
        xge_AddPopup ( popup00 );
        xge_GrabFocus ( popup00, true );
        return true;

  case btn02EDIT:  /* Edit */
        xge_SetMenuWidgets ( menu0, menu00list, true );
        return true;

  case btn02VIEW:  /* View */
        xge_SetMenuWidgets ( menu0, menu01list, true );
        return true;

  case btn02ABOUT:  /* About */
        xge_DisplayInfoMessage ( InfoMsg, -1 );
        return true;

  case btn00REFINE: /* Refine */
        if ( !uniform )
          xge_DisplayErrorMessage ( ErrMsgCannotRefine, -1 );
        else if ( RefineUniform () )
          xge_Redraw ();
        return true;

  case btn00RESET:  /* Reset */
        kwind.panning = cwind.panning = cwind.display_coord =
        kwind.moving_many = kwind.locked = false;
        xge_2DwindInitProjection ( &cwind, -1.0, 1.0, -1.0, 1.0 );
        ResetObject ();
        xge_Redraw ();
        return true;

  case btn00pOPEN:  /* Open */
        xge_RemovePopup ( true );
        SetCurrentDir ();
        xge_SetupFileList ( &filelist, ".", file_filter );
        xge_SetupDirList ( &dirlist, ".", NULL, NULL );
        xge_AddPopup ( popup01 );
        xge_GrabFocus ( popup01, true );
        return true;

  case btn00pSAVE:  /* Save */
        xge_RemovePopup ( true );
process_save_command:
        if ( FilenameCorrect ( filename ) ) {
          if ( !SaveBSCurve ( filename ) )
            xge_DisplayErrorMessage ( ErrMsgCannotSave, -1 );
        }
        else
          goto open_popup02;
        return true;

  case btn00pEXPORT:
        xge_RemovePopup ( true );
        ExportPovRay ();
        return true;

  case btn00pSAVEAS: /* Save as */
        xge_RemovePopup ( true );
open_popup02:
        SetCurrentDir ();
        xge_SetupDirList ( &dirlist, ".", NULL, NULL );
        xge_AddPopup ( popup02 );
        xge_GrabFocus ( popup02, true );
        return true;

  case btn00pEXIT:  /* Exit */
        xge_RemovePopup ( true );
        xge_AddPopup ( popup03 );
        xge_GrabFocus ( popup03, true );
        return true;

  case btn01pOPEN:  /* Open file */
        goto open_the_file;

  case btn01pCANCEL:  /* Cancel opening file */
        xge_RemovePopup ( true );
        xge_ClearListBox ( &filelist );
        xge_ClearListBox ( &dirlist );
        return true;

  case btn02pSAVE:  /* Save file */
        xge_RemovePopup ( true );
        xge_ClearListBox ( &dirlist );
        if ( !SaveBSCurve ( filename ) )
          xge_DisplayErrorMessage ( ErrMsgCannotSave, -1 );
        return true;

  case btn02pCANCEL:  /* Cancel saving file */
        xge_RemovePopup ( true );
        xge_ClearListBox ( &dirlist );
        return true;

  case btn03EXIT:
        xge_done = 1;
        return true;

  case btn03SAVE:
        xge_RemovePopup ( true );
        goto process_save_command;

  case btn03CANCEL:
        xge_RemovePopup ( true );
        return true;

  default:
        break;
      }
      return false;

case xgemsg_SWITCH_COMMAND:
      switch ( er->id ) {
  case sw00SELECT:  /* select/unselect */
        if ( cwind.selecting_mode )
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_SELECTING_TOOL );
        else
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        return true;
  case sw00MOVE:  /* move */
        if ( cwind.moving_tool )
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_MOVING_TOOL );
        else
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        return true;
  case sw00SCALE:  /* scale */
        if ( cwind.scaling_tool )
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_SCALING_TOOL );
        else
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        return true;
  case sw00ROTATE:  /* rotate */
        if ( cwind.rotating_tool )
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_ROTATING_TOOL );
        else
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        return true;
  case sw00SHEAR:  /* shear transformation */
        if ( cwind.shear_tool )
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_SHEAR_TOOL );
        else
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        return true;
  case sw01CURVE:
  case sw01POLYLINE:
  case sw01BEZPOLY:
  case sw01TICKS:
  case sw01CONVH:
  case sw01BASIS:
        goto redraw;
  case sw00NURBS:
        UstawNURBS ();
        funkcja = false;
        xge_Redraw ();
        return true;
  case sw00CLOSED:
        if ( UstawZamknieta () ) {
          RedrawCurveAndKnotsArea ();
          return true;
        }
        else {
          xge_DisplayErrorMessage ( ErrMsgCannotClose, -1 );
          return true;
        }
  case sw00FUNCTION:
        UstawFunkcje ();
        ProjectCurve ();
        nurbs = false;
        xge_Redraw ();
        return false;

  case sw00UNIFORM:  /* Uniform */
        if ( uniform ) {
          SetUniformKnots ();
          xge_Redraw ();
        }
        return true;

  case sw01POLARF:
        UstawPolarf ();
        goto redraw;
  case sw01CURVGR:
        UstawCurvGraph ();
redraw:
        xge_SetClipping ( cwind.er );
        cwind.er->redraw ( cwind.er, true );
        return true;
  case sw00PANZOOM:  /* pan & zoom */
        if ( (kwind.panning = cwind.panning) )
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_PANNING_TOOL );
        else
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        return true;
  case sw00COORD:  /* display coordinates */
        kwind.display_coord = cwind.display_coord;
        return true;
  case sw00MOVEMANY:  /* move many knots */
        return true;
  default:
        break;
      }
      return false;

case xgemsg_SLIDEBAR_COMMAND:
      switch ( er->id ) {
  case sl01CURVGRSC:
        if ( curvgr ) {
          xge_SetClipping ( cwind.er );
          cwind.er->redraw ( cwind.er, true );
        }
        return true;
  default:
        break;
      }
      return false;

case xgemsg_INT_WIDGET_COMMAND:
      switch ( er->id ) {
  case intw00DEG:  /* degree elevation or reduction */
        if ( key > kwind.degree ) {
          if ( key > MAX_DEGREE ) {
            xge_DisplayErrorMessage ( ErrMsgRaiseDeg, -1 );
            return false;
          }
          if ( DegreeElevation () )
            xge_Redraw ();
        }
        else if ( key < kwind.degree ) {
          if ( key < 1 ) {
            xge_DisplayErrorMessage ( ErrMsgReduceDeg, -1 );
            return false;
          }
          if ( DegreeReduction () )
            xge_Redraw ();
        }
        return true;

  case intw01CURVGRDENS:
        curv_graph_dens = key;
        if ( curvgr ) {
          xge_SetClipping ( cwind.er );
          cwind.er->redraw ( cwind.er, true );
        }
        return true;

  default:
        break;
      }
      return false;

case xgemsg_LISTBOX_ITEM_PICK:
      switch ( er->id ) {
  case lb01pDIRLIST:
        if ( !chdir ( &dirlist.itemstr[dirlist.itemind[dirlist.currentitem]] ) ) {
          xge_SetupFileList ( &filelist, ".", file_filter );
          xge_SetupDirList ( &dirlist, ".", NULL, current_directory );
          SetCurrentDir ();
          xge_SetClipping ( popup01 );
          popup01->redraw ( popup01, true );
        }
        return true;

  case lb01pFILELIST:
open_the_file:
        xge_RemovePopup ( true );
        if ( xge_GetCurrentListBoxString ( &filelist, filename ) ) {
          xge_ClearListBox ( &filelist );
          xge_ClearListBox ( &dirlist );
          if ( ReadBSCurve ( filename ) ) {
            ResizeObject ();
/*            ClearPointMarking ();*/
          }
          else {
            ResetObject ();
            xge_DisplayErrorMessage ( ErrMsgCannotRead, -1 );
          }
        }
        xge_Redraw ();
        return true;

  case lb02pDIRLIST:
        if ( !chdir ( &dirlist.itemstr[dirlist.itemind[dirlist.currentitem]] ) ) {
          xge_SetupDirList ( &dirlist, ".", NULL, current_directory );
          SetCurrentDir ();
          xge_SetClipping ( popup02 );
          popup02->redraw ( popup02, true );
        }
        return true;

  default:
        break;
      }
      return false;

case xgemsg_2DWIN_PROJCHANGE:
      ProjectCurve ();
      return true;

case xgemsg_2DWIN_PICK_POINT:
      return FindNearestPoint ( x, y, xge_MINDIST );

case xgemsg_2DWIN_MOVE_POINT:
      SetCPoint ( x, y );
      return true;

case xgemsg_2DWIN_SELECT_POINTS:
      SelectCPoints ( cwind.selection_rect.x0, cwind.selection_rect.x1,
                      cwind.selection_rect.y0, cwind.selection_rect.y1 );
      return true;

case xgemsg_2DWIN_UNSELECT_POINTS:
      UnSelectCPoints ( cwind.selection_rect.x0, cwind.selection_rect.x1,
                        cwind.selection_rect.y0, cwind.selection_rect.y1 );
      return true;

case xgemsg_2DWIN_SAVE_POINTS:
      SaveCPoints ();
      return true;

case xgemsg_2DWIN_TRANSFORM_POINTS:
      TransformCPoints ( &cwind.gwtrans );
      return true;

case xgemsg_2DWIN_FIND_REFBBOX:
      FindBoundingBox ( &cwind.RefBBox );
      return true;

case xgemsg_KNOTWIN_CHANGE_KNOT:
      if ( funkcja )
        UstawFunkcje ();
      if ( uniform ) {
        uniform = false;
        xge_Redraw ();
      }
      else {
        xge_SetClipping ( cwind.er );
        cwind.er->redraw ( cwind.er, true );
      }
      return true;

case xgemsg_KNOTWIN_INSERT_KNOT:
      uu = uniform;
      InsertKnot ();
      if ( uu ) {
        uniform = false;
        xge_Redraw ();
      }
      else {
        xge_SetClipping ( cwind.er );
        cwind.er->redraw ( cwind.er, true );
      }
      return true;

case xgemsg_KNOTWIN_REMOVE_KNOT:
      uu = uniform;
      RemoveKnot ();
      if ( uu ) {
        uniform = false;
        xge_Redraw ();
      }
      else {
        xge_SetClipping ( cwind.er );
        cwind.er->redraw ( cwind.er, true );
      }
      return true;

case xgemsg_KNOTWIN_CHANGE_MAPPING:
      if ( funkcja )
        UstawFunkcje ();
      if ( funkcja || baza ) {
        xge_SetClipping ( cwind.er );
        cwind.er->redraw ( cwind.er, true );
      }
      return true;

case xgemsg_KNOTWIN_ERROR:
      switch ( key ) {
  case 0:
        xge_DisplayErrorMessage ( ErrMsgCannotInsert, -1 );
        return true;
  case 1:
        xge_DisplayErrorMessage ( ErrMsgToManyKnots, -1 );
        return true;
  case 2:
        xge_DisplayErrorMessage ( ErrMsgCannotRemove, -1 );
        return true;
  default:
        break;
      }
      return false;

default:
      return false;
    }
  }
  else {
    switch ( msg ) {
case xgemsg_RESIZE:  /* x and y are new window width and height */
      cwind.er->w = (short)(x-110);      cwind.er->h = (short)(y-78);
      kwind.er->w = cwind.er->w;         kwind.er->y = (short)(y-56);
      menu1->msgproc ( menu1, xgemsg_RESIZE, 0, x, 20 );
      menu0->msgproc ( menu0, xgemsg_RESIZE, 0, 110, (short)(y-20) );
      xge_2DwindInitProjection ( &cwind, cwind.RefBBox.x0, cwind.RefBBox.x1,
                                 cwind.RefBBox.y0, cwind.RefBBox.y1 );
      ResizeObject ();
      xge_Redraw ();
      return true;

case xgemsg_KEY:
      switch ( key ) {
  case 'M':
        xgeResizeWindow ( xge_MAX_WIDTH, xge_MAX_HEIGHT );
        return true;
  case 'm':
        xgeResizeWindow ( xge_WIDTH, xge_HEIGHT );
        return true;
  case 'D': case 'd':
        xgeMoveWindow ( 0, 512 );
        return true;
  case 'S': case 's':
        DumpData ();
        return true;
  default:
        return false;
      }

default:
      break;
    }
  }
  return false;
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

