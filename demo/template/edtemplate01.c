
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "xgedit.h"

#include "spltemplate.h"
#include "edtempwidgets.h"
#include "edtemplate.h"


void init_edwin ( void )
{
  xge_widget *w;

  if ( !pkv_InitScratchMem ( 33554432 ) ) {/* 32MB */
    printf ( "Error: cannot allocate scratch memory stack\n" );
    exit ( 1 );
  }

  /* create widgets for the window 0 */
  win0 = xge_CurrentWindow ();
  xge_New3Dwind ( 0, NULL, win3D0, xge_WIDTH-110, xge_HEIGHT-36, 110, 20,
                  &swind, RysujOkno, RysujPOkno );
    /* the top menu */
  w = xge_NewButton ( 0, NULL, btn00FILE, 76, 19, 0, 0, txtFile );
  w = xge_NewButton ( 0, w, btn00ABOUT, 76, 19, xge_WIDTH-76, 0, txtAbout );
  xge_SetWidgetPositioning ( w, 1, -76, 0 );
  menu0 = xge_NewMenu ( 0, swind.fww.er, MENU0, xge_WIDTH, 20, 0, 0, w);

    /* the side menu */
  w = xge_NewSwitch ( 0, NULL, sw01PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &swind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 0, w, sw01COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &swind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 0, w, sw01aMARK_UNMARK, 109, 16, 0, 21,
                      txtMark_unmark, &swind.selecting_mode );
  w = xge_NewSwitch ( 0, w, sw01aMOVE, 109, 16, 0, 41,
                      txtMove, &swind.moving_tool );
  w = xge_NewSwitch ( 0, w, sw01aSCALE, 109, 16, 0, 61,
                      txtScale, &swind.scaling_tool );
  menu01alist = xge_NewSwitch ( 0, w, sw01aROTATE, 109, 16, 0, 81,
                                txtRotate, &swind.rotating_tool );
  menu1 = xge_NewMenu ( 0, menu0, MENU1, 110, xge_HEIGHT-36, 0, 20,
                        menu01alist );
    /* the bottom menu */
  status0sw = xge_NewSwitch ( 0, NULL, sw02STATUS, 109, 16, 0, xge_HEIGHT-16,
                              txtNULL, &status_0_sw );
  xge_SetWidgetPositioning ( status0sw, 0, 0, 0 );
  status0 = xge_NewTextWidget ( 0, status0sw, STATUSLINE0, xge_MAX_WIDTH-110, 14,
                                110, xge_HEIGHT-14, statustext0 );
  xge_SetWidgetPositioning ( status0, 0, 110, 1 );
  menu2 = xge_NewMenu ( 0, menu1, MENU2, xge_WIDTH, 16, 0, xge_HEIGHT-16, status0 );

    /* the File popup menu */
  w = xge_NewButton ( 0, NULL, btnP00NEW, 76, 19, 2, 22, txtNew );
  w = xge_NewButton ( 0, w, btnP00OPEN, 76, 19, 2, 42, txtOpen );
  w = xge_NewButton ( 0, w, btnP00SAVE, 76, 19, 2, 62, txtSave );
  w = xge_NewButton ( 0, w, btnP00SAVEAS, 76, 19, 2, 82, txtSave_as );
  w = xge_NewButton ( 0, w, btnP00EXIT, 76, 19, 2, 102, txtExit );
  popup00 = xge_NewFMenu ( 0, NULL, POPUP00, 80, 103, 0, 20, w );
  popup00->msgproc = xge_PopupMenuMsg;

    /* the Open file popup menu */
  w = xge_NewTextWidget ( 0, NULL, txtP01DIRNAME, 380, 16, 20+10, 40+10,
                          current_directory );
  w = xge_NewButton ( 0, w, btnP01OPEN, 76, 19, 20+82, 40+180-30, txtOpen );
  w = xge_NewButton ( 0, w, btnP01CANCEL, 76, 19, 20+242, 40+180-30, txtCancel );
  w = xge_NewListBox ( 0, w, lbP01DIRLIST, 180, 99, 20+10, 40+40, &dirlist );
  w = xge_NewListBox ( 0, w, lbP01FILELIST, 180, 99, 220+10, 40+40, &filelist );
  popup01 = xge_NewFMenu ( 0, NULL, POPUP01, 400, 180, 20, 40, w );

    /* the Save as file popup menu */
  w = xge_NewTextWidget ( 0, NULL, txtP02DIRNAME, 380, 16, 20+10, 40+10,
                          current_directory );
  w = xge_NewButton ( 0, w, btnP02SAVE, 76, 19, 20+82, 40+180-30, txtSave );
  w = xge_NewButton ( 0, w, btnP02CANCEL, 76, 19, 20+242, 40+180-30, txtCancel );
  w = xge_NewListBox ( 0, w, lbP02DIRLIST, 180, 99, 20+10, 40+40, &dirlist );
  w = xge_NewTextWidget ( 0, w, txtP02SAVE_AS, 120, 16, 220+10, 40+38, txtSave_as );
  w = xge_NewStringEd ( 0, w, txtedP02FILENAME, 180, 19, 220+10, 40+56,
                        MAX_FILENAME_LGT, filename, &filename_editor );
  popup02 = xge_NewFMenu ( 0, NULL, POPUP02, 400, 180, 20, 40, w );

  xge_SetWinEdRect ( menu2 );
  ResizeWinStatus ( win0 );


  /* create widgets for the window 1 */
  win1 = xge_NewWindow ( "" );
  xge_SetWindow ( win1 );

  xge_New2Dwind ( 1, NULL, win2D0, xge_WIDTH-110, xge_HEIGHT-80, 110, 20,
                  &domwind, RysujDomOkno );
  knwind = xge_NewWidget ( 1, domwind.er, win1D0, xge_WIDTH-110, 40, 110,
                      xge_HEIGHT-56, NULL, NULL, KnotOknoMsg, RysujKnotOkno );
    /* the top menu */
  w = NULL;
  menu3 = xge_NewMenu ( 1, knwind, MENU3, xge_WIDTH, 20, 0, 0, w);

    /* the side menu */
  w = xge_NewSwitch ( 1, NULL, sw11PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &domwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 1, w, sw11COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &domwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 1, w, sw11aMARK_UNMARK, 109, 16, 0, 21,
                      txtMark_unmark, &domwind.selecting_mode );
  w = xge_NewSwitch ( 1, w, sw11aMOVE, 109, 16, 0, 41,
                      txtMove, &domwind.moving_tool );
  w = xge_NewSwitch ( 1, w, sw11aSCALE, 109, 16, 0, 61,
                      txtScale, &domwind.scaling_tool );
  menu14alist = xge_NewSwitch ( 1, w, sw11aROTATE, 109, 16, 0, 81,
                                txtRotate, &domwind.rotating_tool );
  menu4 = xge_NewMenu ( 1, menu3, MENU4, 110, xge_HEIGHT-36, 0, 20,
                        menu14alist );
    /* the bottom menu */
  status1sw = xge_NewSwitch ( 1, NULL, sw15STATUS, 109, 16, 0, xge_HEIGHT-16,
                              txtNULL, &status_1_sw );
  xge_SetWidgetPositioning ( status1sw, 0, 0, 0 );
  status1 = xge_NewTextWidget ( 1, status1sw, STATUSLINE1, xge_MAX_WIDTH-110, 14,
                                110, xge_HEIGHT-14, statustext1 );
  xge_SetWidgetPositioning ( status1, 0, 110, 1 );
  menu5 = xge_NewMenu ( 1, menu4, MENU5, xge_WIDTH, 16, 0, xge_HEIGHT-16, status1 );
  xge_SetWinEdRect ( menu5 );
  ResizeWinStatus ( win1 );


  xge_RedrawAll ();
  xge_SetWindow ( win0 );
  xge_DisplayInfoMessage ( InfoMsg, -1 );
} /*init_edwin*/

void destroy_edwin ( void )
{
  pkv_DestroyScratchMem ();
} /*destroy_edwin*/

