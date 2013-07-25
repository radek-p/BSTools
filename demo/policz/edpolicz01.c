
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
#include "eg1holed.h"
#include "eg2holed.h"
#include "xgedit.h"

#include "render.h"
#include "edpolicz.h"
#include "edwidgets.h"
#include "splhole.h"

void SetupWindow0Widgets ( void )
{
  xge_widget *w;

        /* create widgets for the window 0 */
  win0 = xge_CurrentWindow ();
  xge_New3Dwind ( 0, NULL, win3D0, xge_WIDTH-110, xge_HEIGHT-36, 110, 20,
                  &swind, RysujOkno, RysujPOkno );
    /* the top menu */
  w = xge_NewButton ( 0, NULL, btn00FILE, 58, 19, 0, 0, txtFile );
  w = xge_NewButton ( 0, w, btn00DATA, 58, 19, 60, 0, txtData );
  w = xge_NewButton ( 0, w, btn00EDIT, 58, 19, 120, 0, txtEdit );
  w = xge_NewButton ( 0, w, btn00VIEW, 58, 19, 180, 0, txtView );
  w = xge_NewButton ( 0, w, btn00PICTURE, 58, 19, 240, 0, txtPicture );
  w = xge_NewButton ( 0, w, btn00ABOUT, 58, 19, xge_WIDTH-58, 0, txtAbout );
  xge_SetWidgetPositioning ( w, 1, -58, 0 );
  menu0 = xge_NewMenu ( 0, swind.fww.er, MENU0, xge_WIDTH, 20, 0, 0, w );

    /* the side menu */
      /* Data */
  w = xge_NewSwitch ( 0, NULL, sw01PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &swind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 0, w, sw01COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &swind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewIntWidget ( 0, w, intw01aHOLE_SIDES, 109, 19, 0, 20,
                         3, 16, &hole_k_intw, txtSides, &hole_k );
  w = xge_NewSlidebard ( 0, w, sl01aPARAM0, 109, 10, 0, 42, &surfcparam[0] );
  w = xge_NewSlidebard ( 0, w, sl01aPARAM1, 109, 10, 0, 56, &surfcparam[1] );
  w = xge_NewSlidebard ( 0, w, sl01aPARAM2, 109, 10, 0, 70, &surfcparam[2] );
  w = xge_NewSlidebard ( 0, w, sl01aPARAM3, 109, 10, 0, 84, &surfcparam[3] );
  w = xge_NewSlidebard ( 0, w, sl01aPARAM4, 109, 10, 0, 98, &surfcparam[4] );
  menu01alist = w;
      /* Edit */
  w = xge_NewSwitch ( 0, NULL, sw01PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &swind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 0, w, sw01COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &swind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 0, w, sw01bMARK_UNMARK, 109, 16, 0, 21,
                      txtMark_unmark, &swind.selecting_mode );
  w = xge_NewSwitch ( 0, w, sw01bMOVE, 109, 16, 0, 41,
                      txtMove, &swind.moving_tool );
  w = xge_NewSwitch ( 0, w, sw01bSCALE, 109, 16, 0, 61,
                      txtScale, &swind.scaling_tool );
  w = xge_NewSwitch ( 0, w, sw01bROTATE, 109, 16, 0, 81,
                      txtRotate, &swind.rotating_tool );
  menu01blist = xge_NewSwitch ( 0, w, sw01bSHEAR, 109, 16, 0, 101,
                                txtShear, &swind.shear_tool );
      /* Constraints */
  w = xge_NewSwitch ( 0, NULL, sw01c1ST_SURF_CONSTR, 54, 16, 0, 20,
                      txtFirst, &constraints1 );
  w = xge_NewSwitch ( 0, w, sw01c2ND_SURF_CONSTR, 54, 16, 55, 20,
                      txtSecond, &constraints2 );
  w = xge_NewSwitch ( 0, w, sw01c1ST_CONSTR, 109, 16, 0, 40,
                      txtFirstType, &constraints1st );
  w = xge_NewSwitch ( 0, w, sw01c2ND_CONSTR, 109, 16, 0, 60,
                      txtSecondType, &constraints2nd );
  w = xge_NewSwitch ( 0, w, sw01cZERO_CONSTR, 109, 16, 0, 80,
                      txtZero, &constraintsZero );
  w = xge_NewButton ( 0, w, btn01cGET_CURRENT, 58, 19, 0, 100,
                      txtCurrent );
  if ( !ConstructConstraintSwitches () )
    exit ( 1 );
  menu01cscrolled = xge_NewMenu ( 0, NULL, menu01cSCROLLED_SW,
                                  120, 240, 0, 0, NULL );
  w = xge_NewScrollWidget ( 0, w, scw01cCONSTR_SWITCHES, 110, 74, 0, 122,
                            &scroll_constr_sw, menu01cscrolled );
  ConfigureConstraintWidgets ( true );
  w = xge_NewSwitch ( 0, w, sw01cMARK_UNMARK, 109, 16, 0, xge_HEIGHT-156,
                      txtMark_unmark, &swind.selecting_mode );
  xge_SetWidgetPositioning ( w, 2, 0, -140 );
  w = xge_NewSwitch ( 0, w, sw01cMOVE, 109, 16, 0, xge_HEIGHT-136,
                      txtMove, &swind.moving_tool );
  xge_SetWidgetPositioning ( w, 2, 0, -120 );
  w = xge_NewSwitch ( 0, w, sw01cSCALE, 109, 16, 0, xge_HEIGHT-116,
                      txtScale, &swind.scaling_tool );
  xge_SetWidgetPositioning ( w, 2, 0, -100 );
  w = xge_NewSwitch ( 0, w, sw01cROTATE, 109, 16, 0, xge_HEIGHT-96,
                      txtRotate, &swind.rotating_tool );
  xge_SetWidgetPositioning ( w, 2, 0, -80 );
  w = xge_NewSwitch ( 0, w, sw01cSHEAR, 109, 16, 0, xge_HEIGHT-76,
                      txtShear, &swind.shear_tool );
  xge_SetWidgetPositioning ( w, 2, 0, -60 );
  w = xge_NewSwitch ( 0, w, sw01PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &swind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 0, w, sw01COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &swind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  menu01clist = w;
      /* View */
  w = xge_NewSwitch ( 0, NULL, sw01PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &swind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 0, w, sw01COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &swind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 0, w, sw01dCPOINTS, 109, 16, 0, 21,
                      txtControl_net, &view_surf_cp );
  w = xge_NewSwitch ( 0, w, sw01dSURFACE, 109, 16, 0, 41,
                      txt_surface, &view_surf_spatches );
  w = xge_NewSwitch ( 0, w, sw01dFIRST, 109, 16, 0, 61,
                      txtFirst, &view_surf_1 );
  w = xge_NewSwitch ( 0, w, sw01dSECOND, 109, 16, 0, 81,
                      txtSecond, &view_surf_2 );
  w = xge_NewSwitch ( 0, w, sw01dNUMBERS, 109, 16, 0, 101,
                      txtNumbers, &view_surf_numbers );
  w = xge_NewSwitch ( 0, w, sw01dCONSTRAINT_FRAME, 109, 16, 0, 121,
                      txtConstrFrame, &view_constraints_frame );
  menu01dlist = w;
      /* Picture */
  w = xge_NewSwitch ( 0, NULL, sw01PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &swind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 0, w, sw01COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &swind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = renderbtn0 = xge_NewButton ( 0, w, btn01eRENDER_STOP,
                                   58, 19, 0, 21, txtRender );
  w = xge_NewTextWidget ( 0, w, 0, 12, 16, 0, 43, txtC );
  w = xge_NewTextWidget ( 0, w, 0, 89, 16, 20, 43, txtD_shape_func );
  w = xge_NewSwitch ( 0, w, sw01eGAUSSIAN_C, 16, 16, 0, 61, NULL,
                      &swGaussian_c );
  w = xge_NewSwitch ( 0, w, sw01eMEAN_C, 16, 16, 0, 81, NULL, &swMean_c );
  w = xge_NewSwitch ( 0, w, sw01eLAMBERTISO_C, 16, 16, 0, 101, NULL,
                      &swLambert_c );
  w = xge_NewSwitch ( 0, w, sw01eREFLECTION_C, 16, 16, 0, 121, NULL,
                      &swReflection_c );
  w = xge_NewSwitch ( 0, w, sw01eHIGHLIGHT_C, 16, 16, 0, 141, NULL,
                      &swHighlight_c );
  w = xge_NewSwitch ( 0, w, sw01eSECTIONS_C, 16, 16, 0, 161, NULL,
                      &swSections_c );
  w = xge_NewSwitch ( 0, w, sw01eGAUSSIAN_D, 89, 16, 20, 61,
                      txtGaussian, &swGaussian_d );
  w = xge_NewSwitch ( 0, w, sw01eMEAN_D, 89, 16, 20, 81, txtMean, &swMean_d );
  w = xge_NewSwitch ( 0, w, sw01eLAMBERTISO_D, 89, 16, 20, 101,
                      txtIsophotes, &swLambert_d );
  w = xge_NewSwitch ( 0, w, sw01eREFLECTION_D, 89, 16, 20, 121,
                      txtReflection, &swReflection_d );
  w = xge_NewSwitch ( 0, w, sw01eHIGHLIGHT_D, 89, 16, 20, 141,
                      txtHighlight, &swHighlight_d );
  w = xge_NewSwitch ( 0, w, sw01eSECTIONS_D, 89, 16, 20, 161,
                      txtSections, &swSections_d );
  w = xge_NewSwitch ( 0, w, sw01eSHADOWS, 89, 16, 0, 185,
                      txtShadows, &swShadows );
  w = xge_NewSwitch ( 0, w, sw01eANTIALIAS, 89, 16, 0, 205,
                      txtAntialias, &swAntialias );
  menu01elist = w;

  menu1 = xge_NewMenu ( 0, menu0, MENU1, 110, xge_HEIGHT-36, 0, 20,
                        menu01blist );

    /* Light */
  w = xge_NewSwitch ( 0, NULL, sw01PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &swind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 0, w, sw01COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &swind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = renderbtn1 = xge_NewButton ( 0, w, btn01fRENDER_STOP, 58, 19, 0, 21,
                                   txtRender );
  w = xge_NewSwitch ( 0, w, sw01fLIGHT0DIR, 89, 16, 0, 42, txtLight_0,
                      &edit_light_sw[0] );
  w = xge_NewSlidebard ( 0, w, sl01fLIGHT0INT, 109, 10, 0, 60, &light_int[0] );
  w = xge_NewSwitch ( 0, w, sw01fLIGHT1DIR, 89, 16, 0, 72, txtLight_1,
                      &edit_light_sw[1] );
  w = xge_NewSlidebard ( 0, w, sl01fLIGHT1INT, 109, 10, 0, 90, &light_int[1] );
  w = xge_NewSwitch ( 0, w, sw01fLIGHT2DIR, 89, 16, 0, 102, txtLight_2,
                      &edit_light_sw[2] );
  w = xge_NewSlidebard ( 0, w, sl01fLIGHT1INT, 109, 10, 0, 120, &light_int[2] );
  w = xge_NewSwitch ( 0, w, sw01fLIGHT3DIR, 89, 16, 0, 132, txtLight_3,
                      &edit_light_sw[3] );
  w = xge_NewSlidebard ( 0, w, sl01fLIGHT1INT, 109, 10, 0, 150, &light_int[3] );
  w = xge_NewTextWidget ( 0, w, 0, 109, 16, 0, 162, txtAmbient );
  w = xge_NewSlidebard ( 0, w, sl01fLIGHTAMB, 109, 10, 0, 180, &light_int[4] );
  w = xge_NewSwitch ( 0, w, sw01fREFLECTIONFRAME, 109, 16, 0, 193,
                      txtReflectionFrame, &edit_reflection_frame );
  w = xge_NewSwitch ( 0, w, sw01fHIGHLIGHTFRAME, 109, 16, 0, 213,
                      txtHighlightFrame, &edit_highlight_frame );
  w = xge_NewSwitch ( 0, w, sw01fSECTIONSFRAME, 109, 16, 0, 233,
                      txtSectionsFrame, &edit_sections_frame );
  menu01flist = w;

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

    /* the Light popup menu */
  w = xge_NewButton ( 0, NULL, btnP03SURFACE, 70, 19, 122, 22, txtSurface );
  w = xge_NewButton ( 0, w, btnP03CONSTRAINTS, 70, 19, 122, 42, txtConstraints );
  w = xge_NewButton ( 0, w, btnP03LIGHT, 70, 19, 122, 62, txtLight );
  popup03 = xge_NewFMenu ( 0, NULL, POPUP03, 74, 63, 120, 20, w );
  popup03->msgproc = xge_PopupMenuMsg;

  /* the Exit menu */
  w = xge_NewTextWidget ( 0, NULL, 0, 240, 16, 110, 70, MsgReconsider  );
  w = xge_NewButton ( 0, w, btnP04EXIT, 58, 19, 110, 100, txtExit );
  w = xge_NewButton ( 0, w, btnP04SAVE, 58, 19, 210, 100, txtSave );
  w = xge_NewButton ( 0, w, btnP04CANCEL, 58, 19, 310, 100, txtCancel );
  popup04 = xge_NewFMenu ( 0, NULL, POPUP04, 300, 70, 90, 60, w );
  popup04->msgproc = xge_PopupMenuMsg;

  xge_SetWinEdRect ( menu2 );
  ResizeWinStatus ( win0 );
} /*SetupWindow0Widgets*/

void SetupWindow1Widgets ( void )
{
  xge_widget *w;

        /* create widgets for the window 1 */
  win1 = xge_NewWindow ( "" );
  xge_SetWindow ( win1 );

  xge_New2Dwind ( 1, NULL, win2D0, xge_WIDTH-110, xge_HEIGHT-79, 110, 20,
                  &domwind, RysujDomOkno );
  knotwind = NewGHKnotWind ( 1, domwind.er, win1D0, xge_WIDTH-110, 39, 110,
                      xge_HEIGHT-55, hole_k, &knwind, knots );
    /* the top menu */
  w = xge_NewButton ( 1, NULL, btn10OPTIONS, 58, 19, 0, 0, txtOptions );
  w = xge_NewButton ( 1, w, btn10DATA, 58, 19, 60, 0, txtData );
  w = xge_NewButton ( 1, w, btn10EDIT, 58, 19, 120, 0, txtEdit );
  w = xge_NewButton ( 1, w, btn10VIEW, 58, 19, 180, 0, txtView );
  w = xge_NewButton ( 1, w, btn10INFO, 58, 19, xge_WIDTH-58, 0, txtInfo );
  xge_SetWidgetPositioning ( w, 1, -58, 0 );
  menu3 = xge_NewMenu ( 1, knotwind, MENU3, xge_WIDTH, 20, 0, 0, w );

    /* the side menu */
      /* Options - two menus for two sets of options */
  w = xge_NewSwitch ( 1, NULL, sw11a1FIRST, 54, 16, 0, 21, txtFirst, &sw_opt_1 );
  w = xge_NewSwitch ( 1, w, sw11a1SECOND, 154, 16, 55, 21, txtSecond, &sw_opt_2 );
  w = xge_NewIntWidget ( 1, w, intw11a1ORDER, 109, 19, 0, 41,
                         1, 2, &opt_order_1, txtGCOrder, &options1.order );
  w = xge_NewSwitch ( 1, w, sw11a1RESTRICTED, 109, 16, 0, 63, txtRestricted,
                      &options1.restricted );
  w = xge_NewSwitch ( 1, w, sw11a1COONS, 109, 16, 0, 83, txtCoons,
                      &options1.coons );
  w = xge_NewSwitch ( 1, w, sw11a1BEZIER, 109, 16, 0, 103, txtBezier,
                      &options1.bezier );
  w = xge_NewSwitch ( 1, w, sw11a1SPLINE, 109, 16, 0, 123, txtSpline,
                      &options1.spline );
  w = xge_NewIntWidget ( 1, w, intw11a1NK, 109, 19, 0, 143,
                         1, G2H_S_MAX_NK, &opt_nk_1, txtNknots, &options1.nk );
  w = xge_NewIntWidget ( 1, w, intw11a1M1, 109, 19, 0, 163,
                         1, G2H_S_MAX_M1, &opt_m1_1, txtMult1, &options1.m1 );
  w = xge_NewIntWidget ( 1, w, intw11a1M2, 109, 19, 0, 183,
                         1, G2H_S_MAX_M2, &opt_m2_1, txtMult2, &options1.m2 );
  w = xge_NewSwitch ( 1, w, sw11a1LINEAR, 109, 16, 0, 206, txtLinear,
                      &options1.lin );
  w = xge_NewSwitch ( 1, w, sw11a1QUASIG2, 109, 16, 0, 226, txtQuasiG2,
                      &options1.quasiG2 );
  w = xge_NewSlidebard ( 1, w, sl11a1QUASIG2CONST, 109, 10, 0, 246,
                         &options1.slp );
  w = xge_NewSwitch ( 1, w, sw11a1ALTCENTRE, 109, 16, 0, 263, txtAltCentre,
                      &options1.altcentre );
  w = xge_NewSwitch ( 1, w, sw11a1GAUSSLEGENDRE, 109, 16, 0, 283, txtGaussLegendre,
                      &options1.gausslegendre );
  menu14alist = menu14a1list = w;

  w = xge_NewSwitch ( 1, NULL, sw11a2FIRST, 54, 16, 0, 21, txtFirst, &sw_opt_1 );
  w = xge_NewSwitch ( 1, w, sw11a2SECOND, 54, 16, 55, 21, txtSecond, &sw_opt_2 );
  w = xge_NewIntWidget ( 1, w, intw11a2ORDER, 109, 19, 0, 41,
                         1, 2, &opt_order_2, txtGCOrder, &options2.order );
  w = xge_NewSwitch ( 1, w, sw11a2RESTRICTED, 109, 16, 0, 63, txtRestricted,
                      &options2.restricted );
  w = xge_NewSwitch ( 1, w, sw11a2COONS, 109, 16, 0, 83, txtCoons,
                      &options2.coons );
  w = xge_NewSwitch ( 1, w, sw11a2BEZIER, 109, 16, 0, 103, txtBezier,
                      &options2.bezier );
  w = xge_NewSwitch ( 1, w, sw11a2SPLINE, 109, 16, 0, 123, txtSpline,
                      &options2.spline );
  w = xge_NewIntWidget ( 1, w, intw11a2NK, 109, 19, 0, 143,
                         1, G2H_S_MAX_NK, &opt_nk_2, txtNknots, &options2.nk );
  w = xge_NewIntWidget ( 1, w, intw11a2M1, 109, 19, 0, 163,
                         1, G2H_S_MAX_M1, &opt_m1_2, txtMult1, &options2.m1 );
  w = xge_NewIntWidget ( 1, w, intw11a2M2, 109, 19, 0, 183,
                         1, G2H_S_MAX_M2, &opt_m2_2, txtMult2, &options2.m2 );
  w = xge_NewSwitch ( 1, w, sw11a2LINEAR, 109, 16, 0, 206, txtLinear,
                      &options2.lin );
  w = xge_NewSwitch ( 1, w, sw11a2QUASIG2, 109, 16, 0, 226, txtQuasiG2,
                      &options2.quasiG2 );
  w = xge_NewSlidebard ( 1, w, sl11a2QUASIG2CONST, 109, 10, 0, 246,
                         &options2.slp );
  w = xge_NewSwitch ( 1, w, sw11a2ALTCENTRE, 109, 16, 0, 263, txtAltCentre,
                      &options2.altcentre );
  w = xge_NewSwitch ( 1, w, sw11a2GAUSSLEGENDRE, 109, 16, 0, 283, txtGaussLegendre,
                      &options2.gausslegendre );
  menu14a2list = w;
      /* Data */
  w = xge_NewSwitch ( 1, NULL, sw11PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &domwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 1, w, sw11COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &domwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewIntWidget ( 1, w, intw11bHOLE_SIDES, 109, 19, 0, 20,
                         3, 16, &hole_k_intw, txtSides, &hole_k );
  w = xge_NewSlidebard ( 1, w, sl11bPARAM0, 109, 10, 0, 42, &domcparam[0] );
  w = xge_NewSlidebard ( 1, w, sl11bPARAM1, 109, 10, 0, 56, &domcparam[1] );
  w = xge_NewSlidebard ( 1, w, sl11bPARAM2, 109, 10, 0, 70, &domcparam[2] );
  w = xge_NewSlidebard ( 1, w, sl11bPARAM3, 109, 10, 0, 84, &domcparam[3] );
  menu14blist = w;
      /* Edit */
  w = xge_NewSwitch ( 1, NULL, sw11PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &domwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 1, w, sw11COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &domwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 1, w, sw11cMARK_UNMARK, 109, 16, 0, 21,
                      txtMark_unmark, &domwind.selecting_mode );
  w = xge_NewSwitch ( 1, w, sw11cMOVE, 109, 16, 0, 41,
                      txtMove, &domwind.moving_tool );
  w = xge_NewSwitch ( 1, w, sw11cSCALE, 109, 16, 0, 61,
                      txtScale, &domwind.scaling_tool );
  w = xge_NewSwitch ( 1, w, sw11cROTATE, 109, 16, 0, 81,
                      txtRotate, &domwind.rotating_tool );
  menu14clist = xge_NewSwitch ( 1, w, sw11cSHEAR, 109, 16, 0, 101,
                                txtShear, &domwind.shear_tool );
      /* View */
  w = xge_NewSwitch ( 1, NULL, sw11PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &domwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 1, w, sw11COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &domwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 1, w, sw11dDOMCPOINTS, 109, 16, 0, 21,
                      txtControl_net, &view_dom_cp );
  w = xge_NewSwitch ( 1, w, sw11dDOMSURRPATCHES, 109, 16, 0, 41,
                      txtSurr_patches, &view_dom_spatches );
  w = xge_NewSwitch ( 1, w, sw11dDOMPATCHES1, 109, 16, 0, 61,
                      txtDomain_patches1, &view_dom_patches1 );
  w = xge_NewSwitch ( 1, w, sw11dDOMPATCHES2, 109, 16, 0, 81,
                      txtDomain_patches2, &view_dom_patches2 );
  w = xge_NewSwitch ( 1, w, sw11dDOMNUMBERS, 109, 16, 0, 101,
                      txtNumbers, &view_dom_numbers );
  menu14dlist = w;

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
} /*SetupWindow1Widgets*/

void init_edwin ( void )
{
  if ( !pkv_InitScratchMem ( 134217728 ) ) {/* 128MB */
    printf ( "Error: cannot allocate scratch memory stack\n" );
    exit ( 1 );
  }
  SetupWindow0Widgets ();
  SetupWindow1Widgets ();

    /* setup the geometric data to be processed */
  InitGHObject ( 5 );
  RendInit ();

  xge_RedrawAll ();
  xge_SetWindow ( win0 );
  xge_DisplayInfoMessage ( InfoMsg, -1 );
} /*init_edwin*/

void destroy_edwin ( void )
{
  RendDestroy ();
  if ( domain1 )
    gh_DestroyDomaind ( domain1 );
  if ( domain2 )
    gh_DestroyDomaind ( domain2 );
  pkv_DestroyScratchMem ();
} /*destroy_edwin*/

