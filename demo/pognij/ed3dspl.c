
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
#include <X11/cursorfont.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "xgedit.h"
#include "bsfile.h"

#include "spline3d.h"
#include "ed3dwidgets.h"
#include "ed3dspl.h"

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
      bsf_WriteBSplineCurved ( 3, 4, true, kwind.degree, kwind.lastknot,
                               knots, kwind.closed, &cpoints[0].x,
                               cpmark, NULL );
    }
    else {
      buf = pkv_GetScratchMem ( (kwind.lastknot-kwind.degree)*sizeof(point3d) );
      if ( !buf )
        goto try_nurbs;
      pkv_Selectd ( kwind.lastknot-kwind.degree, 3, 4, 3,
                    &cpoints[0].x, buf );
      bsf_WriteBSplineCurved ( 3, 3, false, kwind.degree, kwind.lastknot,
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
  if ( !cp || !kn || !mkcp )
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
      memcpy ( cpoints, cp, npoints*sizeof(point4d) );
      memcpy ( cpmark, mkcp, npoints );
      if ( spdimen <= 2 ) {
        for ( i = 0; i < npoints; i++ ) {
          cpoints[i].z = 0.0;
          cpoints[i].w = 1.0;
        }
      }
      for ( i = 0; i < npoints; i++ )
        if ( !cpoints[i].w )
          goto failure;
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
void RysujOkno ( xge_widget *er, boolean onscreen )
{
  int  id;

  id = er->id & 0x03;
  xge_DrawGeomWinBackground ( er );
  xge_3DwindDrawCursorPos ( &cwind, id, xge_xx, xge_yy );
  if ( convh )
    DisplayConvh ( id );
  if ( lamana )
    DisplayCPolygon ( id );
  if ( bezpoly )
    DisplayBezPoly ( id );
  if ( lamana )
    DisplayCPoints ( id );
  if ( polarf )
    DisplayPolarf ( id );
  if ( lamana && nurbs )
    DisplayAuxPoints ( id );
  if ( curve )
    DisplayCurve ( id );
  xge_DrawGeomWinSelectionRect ( er, &cwind.selection_rect );
  xge_3DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujOkno*/

void RysujPOkno ( xge_widget *er, boolean onscreen )
{
  int id;

  id = er->id & 0x03;
  xge_DrawGeomWinBackground ( er );
  if ( convh )
    DisplayConvh ( id );
  if ( lamana )
    DisplayCPolygon ( id );
  if ( bezpoly )
    DisplayBezPoly ( id );
  if ( lamana )
    DisplayCPoints ( id );
  if ( polarf )
    DisplayPolarf ( id );
  if ( lamana && nurbs )
    DisplayAuxPoints ( id );
  if ( curvgr || torsgr )
    DisplayCurvGraph ( id );
  if ( curve )
    DisplayCurve ( id );
  xge_DrawGeomWinSelectionRect ( er, &cwind.selection_rect );
  xge_3DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujPOkno*/

/* ///////////////////////////////////////////////////////////////////////// */
void init_edwin ( void )
{
  xge_widget *w;

  pkv_InitScratchMem ( SCRATCH_MEM_SIZE );

    /* setup the four curve windows */
  cwind.fww.er = xge_New3Dwind ( 0, NULL, CURVEWIN,
                           xge_WIDTH-110, xge_HEIGHT-78, 110, 20,
                           &cwind, RysujOkno, RysujPOkno );

    /* setup the knot window */
  kwind.er = xge_NewKnotWind ( 0, cwind.fww.er, KNOTWIN,
                            xge_WIDTH-110, 56, 110, xge_HEIGHT-56,
                            &kwind, MAX_KNOTS, knots );

    /* setup the top menu */
  w = xge_NewButton ( 0, NULL, btn00FILE, 60, 19, 0, 0, txtFile );
  w = xge_NewButton ( 0, w, btn00EDIT, 60, 19,   62, 0, txtEdit );
  w = xge_NewButton ( 0, w, btn00VIEW, 60, 19,  124, 0, txtView );
  w = xge_NewButton ( 0, w, btn00ABOUT, 60, 19, xge_WIDTH-60, 0, txtAbout );
  xge_SetWidgetPositioning ( w, 1, -60, 0 );
  menu1 = xge_NewMenu ( 0, kwind.er, 5, xge_WIDTH, 20, 0, 0, w );

    /* setup the side menu 00 widgets */
  w = xge_NewIntWidget ( 0, NULL, intw01aDEGREE, 109, 19, 0, 20, 0, MAX_DEGREE+1,
                         &ddegree, txtDegree, &kwind.degree );
  w = xge_NewSwitch ( 0, w, sw01aSELECT, 109, 16, 0,  40, txtSelectUnselect, &cwind.selecting_mode );
  w = xge_NewSwitch ( 0, w, sw01aMOVE, 109, 16, 0,  60, txtMove, &cwind.moving_tool );
  w = xge_NewSwitch ( 0, w, sw01aSCALE, 109, 16, 0,  80, txtScale, &cwind.scaling_tool );
  w = xge_NewSwitch ( 0, w, sw01aROTATE, 109, 16, 0, 100, txtRotate, &cwind.rotating_tool );
  w = xge_NewSwitch ( 0, w, sw01aSHEAR, 109, 16, 0, 120, txtShear, &cwind.shear_tool );
  w = xge_NewSwitch ( 0, w, sw01aNURBS, 109, 16, 0, 140, txtNURBS, &nurbs );
  w = xge_NewSwitch ( 0, w, sw01aCLOSED, 109, 16, 0, 160, txtClosed, &kwind.closed );
  w = xge_NewSwitch ( 0, w, sw01aUNIFORM, 60, 16, 0, 180, txtUniform, &uniform );
  w = xge_NewButton ( 0, w, btn01aREFINE, 60, 19, 0, 200, txtRefine );
  w = xge_NewButton ( 0, w, btn01aRESET, 60, 19, 0, 220, txtReset );
  w = xge_NewSwitch ( 0, w, sw01aPANZOOM, 109, 16, 0, xge_HEIGHT-56,
                           txtPanZoom, &cwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( 0, w, sw01aCOORD, 109, 16, 0, xge_HEIGHT-36,
                           txtCoordinates, &cwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  menu00list = w = xge_NewSwitch ( 0, w, sw01aMOVEMANY, 109, 16, 0, xge_HEIGHT-16,
                                   txtMoveManyKnots, &kwind.moving_many );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );

    /* setup the side menu 01 widgets */
  w = xge_NewSwitch ( 0, NULL, sw01bCPOLY, 109, 16, 0, 20, txtControlPolygon, &lamana );
  w = xge_NewSwitch ( 0, w, sw01bCURVE, 109, 16, 0, 40, txtCurve, &curve );
  w = xge_NewSwitch ( 0, w, sw01bTICKS, 109, 16, 0, 60, txtTicks, &ticks );
  w = xge_NewSwitch ( 0, w, sw01bBEZPOLY, 109, 16, 0, 80, txtBezierPolygons, &bezpoly );
  w = xge_NewSwitch ( 0, w, sw01bCONVH, 109, 16, 0, 100, txtConvexHulls, &convh );
  w = xge_NewSwitch ( 0, w, sw01bDIAGF, 109, 16, 0, 120, txtDiagonalForms, &polarf );
  w = xge_NewSwitch ( 0, w, sw01bCURVGRAPH, 109, 16, 0, 140, txtCurvatureGraph, &curvgr );
  w = xge_NewSlidebard ( 0, w, sl01bCURVSCALE, 109, 10, 0, 160, &curvgraphscale );
  w = xge_NewSwitch ( 0, w, sw01bTORSGRAPH, 109, 16, 0, 174, txtTorsionGraph, &torsgr );
  w = xge_NewSlidebard ( 0, w, sl01bTORSSCALE, 109, 10, 0, 194, &torsgraphscale );
  w = xge_NewIntWidget ( 0, w, intw01bGRAPHDENS, 109, 19, 0, 208, 1, 32,
                         &ccdens, txtGraphDensity, &curv_graph_dens );
  w = xge_NewSwitch ( 0, w, sw01aPANZOOM, 109, 16, 0, xge_HEIGHT-56,
                           txtPanZoom, &cwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( 0, w, sw01aCOORD, 109, 16, 0, xge_HEIGHT-36,
                           txtCoordinates, &cwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  menu01list = w = xge_NewSwitch ( 0, w, sw01aMOVEMANY, 109, 16, 0, xge_HEIGHT-16,
                                   txtMoveManyKnots, &kwind.moving_many );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );

  menu0 = xge_NewMenu ( 0, menu1, 6, 110, xge_HEIGHT-20, 0, 20, menu00list );

  w = xge_NewButton ( 0, NULL, btn00pOPEN, 76, 19, 2, 22, txtOpen );
  w = xge_NewButton ( 0, w, btn00pSAVE, 76, 19, 2, 42, txtSave );
  w = xge_NewButton ( 0, w, btn00pSAVEAS, 76, 19, 2, 62, txtSaveAs );
  w = xge_NewButton ( 0, w, btn00pEXIT, 76, 19, 2, 82, txtExit );
  popup00 = xge_NewFMenu ( 0, NULL, POPUP00, 80, 83, 0, 20, w );
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

/* ////////////////////////////////////////////////////////////////////////// */
int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  int i;
  boolean uu;

  if ( er ) {
    switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
      switch ( er->id ) {
  case btn00FILE:   /* File */
        xge_AddPopup ( popup00 );
        xge_GrabFocus ( popup00, true );
        return true;

  case btn00pEXIT:  /* Exit */
        xge_RemovePopup ( true );   
        xge_AddPopup ( popup03 );
        xge_GrabFocus ( popup03, true );
        return true;

  case btn00EDIT:  /* Edit */
        xge_SetMenuWidgets ( menu0, menu00list, true );
        return true;

  case btn00VIEW:  /* View */
        xge_SetMenuWidgets ( menu0, menu01list, true );
        return true;

  case btn00ABOUT:  /* About */
        xge_DisplayInfoMessage ( InfoMsg, -1 );
        return true;

  case btn01aREFINE: /* Refine */
        if ( !uniform )
          xge_DisplayErrorMessage ( ErrMsgCannotRefine, -1 );
        else if ( RefineUniform () )
          xge_Redraw ();
        return true;

  case btn01aRESET:  /* Reset */
        ResetObject ();
        xge_Redraw ();
        return true;

  case btn00pOPEN:  /* Open */
        xge_RemovePopup ( true );
        SetCurrentDir ();
        xge_SetupFileList ( &filelist, ".", file_filter, false );
        xge_SetupDirList ( &dirlist, ".", NULL, false, NULL );
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

  case btn00pSAVEAS: /* Save as */
        xge_RemovePopup ( true );
open_popup02:
        SetCurrentDir ();
        xge_SetupDirList ( &dirlist, ".", NULL, false, NULL );
        xge_AddPopup ( popup02 );
        xge_GrabFocus ( popup02, true );
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
      break;

case xgemsg_SWITCH_COMMAND:
      switch ( er->id ) {
  case sw01aSELECT:  /* select/unselect */
        if ( cwind.selecting_mode )
          xge_3DwindEnableGeomWidget ( &cwind, xge_3DWIN_SELECTING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &cwind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        return true;
  case sw01aMOVE:  /* move */
        if ( cwind.moving_tool )
          xge_3DwindEnableGeomWidget ( &cwind, xge_3DWIN_MOVING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &cwind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        return true;
  case sw01aSCALE:  /* scale */
        if ( cwind.scaling_tool )
          xge_3DwindEnableGeomWidget ( &cwind, xge_3DWIN_SCALING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &cwind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        return true;
  case sw01aROTATE:  /* rotate */
        if ( cwind.rotating_tool )
          xge_3DwindEnableGeomWidget ( &cwind, xge_3DWIN_ROTATING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &cwind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        return true;
  case sw01aSHEAR:  /* shear */
        if ( cwind.shear_tool )
          xge_3DwindEnableGeomWidget ( &cwind, xge_3DWIN_SHEAR_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &cwind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        return true;
  case sw01aUNIFORM:  /* Uniform */
        if ( uniform ) {
          SetUniformKnots ();
          xge_Redraw ();
        }
        return true;
  case sw01bCPOLY:
  case sw01bCURVE:
  case sw01bTICKS:
  case sw01bBEZPOLY:
  case sw01bCONVH:
        xge_SetClipping ( cwind.fww.er );
        cwind.fww.er->redraw ( cwind.fww.er, true );
        return true;
  case sw01aNURBS:
        UstawNURBS ();
        xge_SetClipping ( cwind.fww.er );
        cwind.fww.er->redraw ( cwind.fww.er, true );
        return true;
  case sw01aCLOSED:
        if ( UstawZamknieta () ) {
          xge_Redraw ();
          return true;
        }
        else
          return false;
  case sw01bDIAGF:
        UstawPolarf ();
        xge_SetClipping ( cwind.fww.er );
        cwind.fww.er->redraw ( cwind.fww.er, true );
        return true;
  case sw01bCURVGRAPH:
        UstawCurvGraph ();
        xge_SetClipping ( cwind.fww.er );
        cwind.fww.er->redraw ( cwind.fww.er, true );
        return true;
  case sw01bTORSGRAPH:
        UstawTorsGraph ();
        xge_SetClipping ( cwind.fww.er );
        cwind.fww.er->redraw ( cwind.fww.er, true );
        return true;
  case sw01aPANZOOM:  /* pan & zoom */
        if ( (kwind.panning = cwind.panning) )
          xge_3DwindEnableGeomWidget ( &cwind, xge_3DWIN_PANNING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &cwind, xge_3DWIN_NO_TOOL );
        xge_Redraw ();
        return true;
  case sw01aCOORD:  /* coordinates */
        kwind.display_coord = cwind.display_coord;
        return true;
  default:
        return false;
      }
      xge_SetClipping ( cwind.fww.er );
      cwind.fww.er->redraw ( cwind.fww.er, true );
      return true;

case xgemsg_SLIDEBAR_COMMAND:
      switch ( er->id ) {
  case sl01bCURVSCALE:
        if ( curvgr ) {
          xge_SetClipping ( cwind.fww.er );
          cwind.fww.er->redraw ( cwind.fww.er, true );
        }
        return true;
  case sl01bTORSSCALE:
        if ( torsgr ) {
          xge_SetClipping ( cwind.fww.er );
          cwind.fww.er->redraw ( cwind.fww.er, true );
        }
        return true;
  default:
        break;
      }
      break;

case xgemsg_INT_WIDGET_COMMAND:
      switch ( er->id ) {
  case intw01aDEGREE: /* degree elevation or reduction */
        uu = uniform;
        if ( key > kwind.degree ) {
          if ( key > MAX_DEGREE ) {
            xge_DisplayErrorMessage ( "Error: Cannot raise degree above the limit.", -1 );
            return false;
          }
          if ( DegreeElevation () ) {
            if ( uu ) {
              uniform = false;
              xge_Redraw ();
            }
            else {
              xge_SetClipping ( kwind.er );
              kwind.er->redraw ( kwind.er, true );
              xge_SetClipping ( cwind.fww.er );
              cwind.fww.er->redraw ( cwind.fww.er, true );
            }
          }
        }
        else if ( key < kwind.degree ) {
          if ( key < 1 ) {
            xge_DisplayErrorMessage ( "Error: Cannot reduce degree below 1.", -1 );
            return false;
          }
          if ( DegreeReduction () ) {
            if ( uu ) {
              uniform = false;
              xge_Redraw ();
            }
            else {
              xge_SetClipping ( kwind.er );
              kwind.er->redraw ( kwind.er, true );
              xge_SetClipping ( cwind.fww.er );
              cwind.fww.er->redraw ( cwind.fww.er, true );
            }
          }
        }
        return true;
  case intw01bGRAPHDENS:
        curv_graph_dens = key;
        if ( curvgr || torsgr ) {
          xge_SetClipping ( cwind.fww.er );
          cwind.fww.er->redraw ( cwind.fww.er, true );
        }
        return true;
  default:
        break;
      }
      break;

case xgemsg_LISTBOX_ITEM_PICK:
      switch ( er->id ) {
  case lb01pDIRLIST:
        if ( !chdir ( &dirlist.itemstr[dirlist.itemind[dirlist.currentitem]] ) ) {
          xge_SetupFileList ( &filelist, ".", file_filter, false );
          xge_SetupDirList ( &dirlist, ".", NULL, false, current_directory );
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
/*            ClearPointMarking (); */
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
          xge_SetupDirList ( &dirlist, ".", NULL, false, current_directory );
          SetCurrentDir ();
          xge_SetClipping ( popup02 );
          popup02->redraw ( popup02, true );
        }
        return true;

  default:
        break;
      }
      break;

case xgemsg_3DWIN_RESIZE:
case xgemsg_3DWIN_PROJCHANGE:
      for ( i = 0; i < 4; i++ )
        ProjectCurve ( i );
      return true;

case xgemsg_3DWIN_PICK_POINT:
      return FindNearestPoint ( er->id, x, y, xge_MINDIST );

case xgemsg_3DWIN_MOVE_POINT:
      SetCPoint ( er->id, x, y );
      return true;

case xgemsg_3DWIN_SELECT_POINTS:
      SelectCPoints ( er->id,
                      cwind.selection_rect.x0, cwind.selection_rect.x1,
                      cwind.selection_rect.y0, cwind.selection_rect.y1 );
      return true;

case xgemsg_3DWIN_UNSELECT_POINTS:
      UnSelectCPoints ( er->id,
                        cwind.selection_rect.x0, cwind.selection_rect.x1,
                        cwind.selection_rect.y0, cwind.selection_rect.y1 );
      return true;

case xgemsg_3DWIN_SAVE_POINTS:
      SaveCPoints ();
      return true;

case xgemsg_3DWIN_TRANSFORM_POINTS:
      TransformCPoints ( &cwind.gwtrans );
      return true;

case xgemsg_3DWIN_FIND_REFBBOX:
      FindBoundingBox ( &cwind.RefBBox );
      return true;

case xgemsg_KNOTWIN_CHANGE_KNOT:
      if ( uniform ) {
        uniform = false;
        xge_Redraw ();
      }
      else {
        xge_SetClipping ( cwind.fww.er );
        cwind.fww.er->redraw ( cwind.fww.er, true );
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
        xge_SetClipping ( cwind.fww.er );
        cwind.fww.er->redraw ( cwind.fww.er, true );
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
        xge_SetClipping ( cwind.fww.er );
        cwind.fww.er->redraw ( cwind.fww.er, true );
      }
      return true;

case xgemsg_KNOTWIN_CHANGE_MAPPING:
      return true;

case xgemsg_KNOTWIN_ERROR:
      switch ( key ) {
  case 0:
        xge_DisplayErrorMessage ( "Error: Cannot insert a knot at this point.", -1 );
        return true;
  case 1:
        xge_DisplayErrorMessage ( "Error: Too many knots.", -1 );
        return true;
  case 2:
        xge_DisplayErrorMessage ( "Error: Cannot remove this knot.", -1 );
        return true;
  default:
        break;
      }
      return false;

default:
      break;
    }
  }
  else {
    switch ( msg ) {
case xgemsg_RESIZE:
      cwind.fww.er->msgproc ( cwind.fww.er, xgemsg_RESIZE, 0,
                (short)(xge_current_width-110), (short)(xge_current_height-78) );
      kwind.er->y = (short)(y-56);
      kwind.er->msgproc ( kwind.er, xgemsg_RESIZE, 0, (short)(x-110), kwind.er->h );
      menu1->msgproc ( menu1, xgemsg_RESIZE, 0, x, 20 );
      menu0->msgproc ( menu0, xgemsg_RESIZE, 0, 110, (short)(y-20) );
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
  default:
        break;
      }
      return false;

default:
      break;
    }
  }
  return false;
} /*CallBack*/

/* ////////////////////////////////////////////////////////////////////////// */
int main ( int argc, char *argv[] )
{
  xge_Init ( argc, argv, CallBack, NULL );
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();
  xge_Cleanup ();
  exit ( 0 );
} /*main*/

