
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <signal.h>

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
#include "camerad.h"
#include "multibs.h"
#include "xgedit.h"

#include "render.h"
#include "spl3d.h"
#include "ed3ds.h"
#include "ed3dswidgets.h"
#include "pomnijipc.h"

/* ///////////////////////////////////////////////////////////////////////// */
void SetupWindow0Widgets ( void )
{
  xge_widget *w;

  /* setup the first window and its widgets */
  win0 = xge_CurrentWindow ();
  w = xge_New3Dwind ( 0, NULL, win3D0, xge_WIDTH-110, xge_HEIGHT-20, 110, 20,
                      &swind, RysujOkno, RysujPOkno );
    /* the top menu */
  w = xge_NewButton ( 0, NULL, btn01FILEpopup, 58, 19, 0, 0, txtFile );
  w = xge_NewButton ( 0, w, btn01EDITpopup, 58, 19, 60, 0, txtEdit );
  w = xge_NewButton ( 0, w, btn01VIEW, 58, 19, 120, 0, txtView );
  w = xge_NewButton ( 0, w, btn01PICTURE, 58, 19, 180, 0, txtPicture );
  w = xge_NewButton ( 0, w, btn01ABOUT, 58, 19, xge_WIDTH-58, 0, txtAbout );
  xge_SetWidgetPositioning ( w, 1, -58, 0 );
  menu1 = xge_NewMenu ( 0, swind.fww.er, MENU1, xge_WIDTH, 20, 0, 0, w );

    /* the side menu 00 - Edit surface */
  w = xge_NewSwitch ( 0, NULL, sw00aPAN_ZOOM, 109, 16, 0, xge_HEIGHT-40,
                                 txtPan_zoom, &swind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 0, w, sw00aCOORDINATES, 109, 16, 0, xge_HEIGHT-20,
                      txtCoordinates, &swind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewButton ( 0, w, btn00aRESET, 58, 19, 0, 21, txtReset );
  w = xge_NewSwitch ( 0, w, sw00aMARK, 109, 16, 0, 41, txtMark_unmark,
                      &swind.selecting_mode );
  w = xge_NewSwitch ( 0, w, sw00aMOVE, 109, 16, 0, 61, txtMove,
                      &swind.moving_tool );
  w = xge_NewSwitch ( 0, w, sw00aSCALE, 109, 16, 0, 81, txtScale,
                      &swind.scaling_tool );
  w = xge_NewSwitch ( 0, w, sw00aROTATE, 109, 16, 0, 101,
                      txtRotate, &swind.rotating_tool );
  menu00list = xge_NewSwitch ( 0, w, sw00aSHEAR, 109, 16, 0, 121,
                               txtShear, &swind.shear_tool );

    /* the side menu 01 - View */
  w = xge_NewSwitch ( 0, NULL, sw00bPAN_ZOOM, 109, 16, 0, xge_HEIGHT-40,
                                 txtPan_zoom, &swind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 0, w, sw00bCOORDINATES, 109, 16, 0, xge_HEIGHT-20,
                      txtCoordinates, &swind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 0, w, sw00bCONTROLNET, 109, 16, 0, 21,
                      txtControl_net, &display_control_net );
  w = xge_NewSwitch ( 0, w, sw00bSURFACE, 109, 16, 0, 41,
                      txtSurface, &display_surface );
  w = xge_NewIntWidget ( 0, w, intw00bU_DENSITY, 80, 19, 0, 60, 1, 16,
                         &ddens_u, txtU_density, &display_bez_dens_u );
  w = xge_NewIntWidget ( 0, w, intw00bV_DENSITY, 80, 19, 0, 81, 1, 16,
                         &ddens_v, txtV_density, &display_bez_dens_v );
  w = xge_NewSwitch ( 0, w, sw00bBEZIER_NETS, 109, 16, 0, 102,
                      txtBezier_nets, &display_Bezier_nets );
  w = xge_NewSwitch ( 0, w, sw00bCONSTRPOLY, 109, 16, 0, 122,
                      txtConstrPoly, &display_constr_poly );
  menu01list = xge_NewSwitch ( 0, w, sw00bTRANSFNET, 109, 16, 0, 142,
                               txtTransfNet, &display_trans_net );

  menu0 = xge_NewMenu ( 0, menu1, MENU0, 110, xge_HEIGHT-56, 0, 20, menu00list );

    /* the side menu 02 - Picture */
  w = xge_NewSwitch ( 0, NULL, sw00aPAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &swind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 0, w, sw00aCOORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &swind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = renderbtn0 = xge_NewButton ( 0, w, btn00cRENDER_STOP,   
                                   58, 19, 0, 21, txtRender );
  w = xge_NewTextWidget ( 0, w, 0, 12, 16, 0, 43, txtC );
  w = xge_NewTextWidget ( 0, w, 0, 89, 16, 20, 43, txtD_shape_func );
  w = xge_NewSwitch ( 0, w, sw00cGAUSSIAN_C, 16, 16, 0, 61, NULL,
                      &swGaussian_c );
  w = xge_NewSwitch ( 0, w, sw00cMEAN_C, 16, 16, 0, 81, NULL, &swMean_c );
  w = xge_NewSwitch ( 0, w, sw00cLAMBERTISO_C, 16, 16, 0, 101, NULL,
                      &swLambert_c );
  w = xge_NewSwitch ( 0, w, sw00cREFLECTION_C, 16, 16, 0, 121, NULL,
                      &swReflection_c );
  w = xge_NewSwitch ( 0, w, sw00cHIGHLIGHT_C, 16, 16, 0, 141, NULL,
                      &swHighlight_c );
  w = xge_NewSwitch ( 0, w, sw00cSECTIONS_C, 16, 16, 0, 161, NULL,
                      &swSections_c );
  w = xge_NewSwitch ( 0, w, sw00cPARAM_C, 16, 16, 0, 181, NULL,
                      &swParam_c );
  w = xge_NewSwitch ( 0, w, sw00cGAUSSIAN_D, 89, 16, 20, 61,
                      txtGaussian, &swGaussian_d );
  w = xge_NewSwitch ( 0, w, sw00cMEAN_D, 89, 16, 20, 81, txtMean, &swMean_d );
  w = xge_NewSwitch ( 0, w, sw00cLAMBERTISO_D, 89, 16, 20, 101,
                      txtIsophotes, &swLambert_d );
  w = xge_NewSwitch ( 0, w, sw00cREFLECTION_D, 89, 16, 20, 121,
                      txtReflection, &swReflection_d );
  w = xge_NewSwitch ( 0, w, sw00cHIGHLIGHT_D, 89, 16, 20, 141,
                      txtHighlight, &swHighlight_d );
  w = xge_NewSwitch ( 0, w, sw00cSECTIONS_D, 89, 16, 20, 161,
                      txtSections, &swSections_d );
  w = xge_NewSwitch ( 0, w, sw00cPARAM_D, 89, 16, 20, 181,
                      txtParamQual, &swParam_d );
  w = xge_NewSlidebar2d ( 0, w, sl00cREND_CFRANGE, 109, 10, 0, 201, cf_range );
  w = xge_NewSlidebard ( 0, w, sl00cREND_DFSF, 109, 10, 0, 215, &dfs_factor );
  w = xge_NewSwitch ( 0, w, sw00cSHADOWS, 89, 16, 0, 229,
                      txtShadows, &swShadows );
  w = xge_NewSwitch ( 0, w, sw00cANTIALIAS, 89, 16, 0, 249,
                      txtAntialias, &swAntialias );
  menu02list = w;

    /* the side menu 03 - Edit light */
  w = xge_NewSwitch ( 0, NULL, sw00aPAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &swind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 0, w, sw00aCOORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &swind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = renderbtn1 = xge_NewButton ( 0, w, btn00dRENDER_STOP, 58, 19, 0, 21,
                                   txtRender );
  w = xge_NewSwitch ( 0, w, sw00dLIGHT0DIR, 89, 16, 0, 42, txtLight_0,
                      &edit_light_sw[0] );
  w = xge_NewSlidebard ( 0, w, sl00dLIGHT0INT, 109, 10, 0, 60, &light_int[0] );
  w = xge_NewSwitch ( 0, w, sw00dLIGHT1DIR, 89, 16, 0, 72, txtLight_1,
                      &edit_light_sw[1] );
  w = xge_NewSlidebard ( 0, w, sl00dLIGHT1INT, 109, 10, 0, 90, &light_int[1] );
  w = xge_NewSwitch ( 0, w, sw00dLIGHT2DIR, 89, 16, 0, 102, txtLight_2,
                      &edit_light_sw[2] );
  w = xge_NewSlidebard ( 0, w, sl00dLIGHT1INT, 109, 10, 0, 120, &light_int[2] );
  w = xge_NewSwitch ( 0, w, sw00dLIGHT3DIR, 89, 16, 0, 132, txtLight_3,
                      &edit_light_sw[3] );
  w = xge_NewSlidebard ( 0, w, sl00dLIGHT1INT, 109, 10, 0, 150, &light_int[3] );
  w = xge_NewTextWidget ( 0, w, 0, 109, 16, 0, 162, txtAmbient );
  w = xge_NewSlidebard ( 0, w, sl00dLIGHTAMB, 109, 10, 0, 180, &light_int[4] );
  w = xge_NewSwitch ( 0, w, sw00dREFLECTIONFRAME, 109, 16, 0, 193,
                      txtReflectionFrame, &edit_reflection_frame );
  w = xge_NewSwitch ( 0, w, sw00dHIGHLIGHTFRAME, 109, 16, 0, 213,  
                      txtHighlightFrame, &edit_highlight_frame );  
  w = xge_NewSwitch ( 0, w, sw00dSECTIONSFRAME, 109, 16, 0, 233,   
                      txtSectionsFrame, &edit_sections_frame );    
  menu03list = w;

  /* setup the File popup menu */
  w = xge_NewButton ( 0, NULL, btn02OPENpopup, 76, 19, 2, 22, txtOpen );
  w = xge_NewButton ( 0, w, btn02SAVEpopup, 76, 19, 2, 42, txtSave );
  w = xge_NewButton ( 0, w, btn02SAVEASpopup, 76, 19, 2, 62, txtSave_as );
  w = xge_NewButton ( 0, w, btn02EXPORT, 76, 19, 2, 82, txtExport );
  w = xge_NewButton ( 0, w, btn02EXIT, 76, 19, 2, 102, txtExit );
  popup00 = xge_NewFMenu ( 0, NULL, POPUP00, 80, 103, 0, 20, w );
  popup00->msgproc = xge_PopupMenuMsg;

  /* setup the Open file menu */
  w = xge_NewTextWidget ( 0, NULL, 0, 380, 16, 20+10, 40+10, current_directory );
  w = xge_NewButton ( 0, w, btn03OPEN, 76, 19, 20+82, 40+180-30, txtOpen );
  w = xge_NewButton ( 0, w, btn03CANCEL, 76, 19, 20+242, 40+180-30, txtCancel );
  w = xge_NewListBox ( 0, w, lb03DIRLIST, 180, 99, 20+10, 40+40, &dirlist );
  w = xge_NewListBox ( 0, w, lb03FILELIST, 180, 99, 220+10, 40+40, &filelist );
  popup01 = xge_NewFMenu ( 0, NULL, POPUP01, 400, 180, 20, 40, w );

  /* setup the Save as file menu */
  w = xge_NewTextWidget ( 0, NULL, 0, 380, 16, 20+10, 40+10, current_directory );
  w = xge_NewButton ( 0, w, btn04SAVE, 76, 19, 20+82, 40+180-30, txtSave );
  w = xge_NewButton ( 0, w, btn04CANCEL, 76, 19, 20+242, 40+180-30, txtCancel );
  w = xge_NewListBox ( 0, w, lb04DIRLIST, 180, 99, 20+10, 40+40, &dirlist );
  w = xge_NewTextWidget ( 0, w, txt04SAVE_AS, 120, 16, 220+10, 40+38, txtSave_as );
  w = xge_NewStringEd ( 0, w, txted04FILENAME, 180, 19, 220+10, 40+56,
                        MAX_FILENAME_LGT, filename, &filename_editor );
  popup02 = xge_NewFMenu ( 0, NULL, POPUP02, 400, 180, 20, 40, w );

  /* setup the Edit menu */
  w = xge_NewButton ( 0, NULL, btn05SURFACE, 76, 19, 62, 22, txt_Surface );
  w = xge_NewButton ( 0, w, btn05LIGHT, 76, 19, 62, 42, txt_Light );
  popup03 = xge_NewFMenu ( 0, NULL, POPUP03, 74, 43, 60, 20, w );
  popup03->msgproc = xge_PopupMenuMsg;

  /* setup the Exit menu */
  w = xge_NewTextWidget ( 0, NULL, 0, 240, 16, 110, 70, MsgReconsider  );
  w = xge_NewButton ( 0, w, btn06EXIT, 58, 19, 110, 100, txtExit );
  w = xge_NewButton ( 0, w, btn06SAVE, 58, 19, 210, 100, txtSave );
  w = xge_NewButton ( 0, w, btn06CANCEL, 58, 19, 310, 100, txtCancel );
  popup04 = xge_NewFMenu ( 0, NULL, POPUP04, 300, 70, 90, 60, w );
  popup04->msgproc = xge_PopupMenuMsg;

    /* the bottom menu */
  status0sw = xge_NewSwitch ( 0, NULL, sw02STATUS, 109, 16, 0, xge_HEIGHT-16,
                              txtNULL, &status_0_sw );
  xge_SetWidgetPositioning ( status0sw, 0, 0, 0 );
  status0 = xge_NewTextWidget ( 0, status0sw, STATUSLINE0, xge_MAX_WIDTH-110, 14,
                                110, xge_HEIGHT-14, statustext0 );
  xge_SetWidgetPositioning ( status0, 0, 110, 1 );
  menu2 = xge_NewMenu ( 0, menu0, MENU2, xge_WIDTH, 16, 0, xge_HEIGHT-16, status0 );
  xge_SetWinEdRect ( menu2 );
  ResizeWinStatus ( win0 );
} /*SetupWindow0Widgets*/

void SetupWindow1Widgets ( void )
{
  xge_widget *w;

  /* setup the second window and its widgets */
  win1 = xge_NewWindow ( "" );
  xge_SetWindow ( win1 );
  domwind = xge_NewT2KnotWind ( 1, NULL, 4, xge_WIDTH-110, xge_HEIGHT-20, 110, 20,
                                50, &kwind, RysujDOkno,
                                MAX_KNOTS, knots_u, MAX_KNOTS, knots_v );
  xge_T2KnotWindSetAltKnots ( &kwind,
            MAX_BLENDING_CONSTRAINTS, 0, 1, blending_constr_knots,
            0, 0, 0, NULL );
    /* the top menu */
  w = xge_NewButton ( 1, NULL, btn14aTYPE, 76, 19, 0, 0, txtSurface_Type );
  menu3 = xge_NewMenu ( 1, domwind, MENU3, xge_WIDTH, 20, 0, 0, w );
    /* the side menu */
  w = xge_NewTextWidget ( 1, NULL, txt13GENERAL, 109, 16, 0, 22,
                          txtGeneral_BSpline );
  w = xge_NewIntWidget ( 1, w, intw13DEGREE_U, 109, 19, 0, 40, 0,
                         MAX_DEGREE+1, &ddeg_u, txtDegree_u, &degree_u );
  w = xge_NewIntWidget ( 1, w, intw13DEGREE_V, 109, 19, 0, 61, 0,
                         MAX_DEGREE+1, &ddeg_v, txtDegree_v, &degree_v );
  w = xge_NewSwitch ( 1, w, sw13CLOSED_U, 109, 16, 0, 82,
                      txtClosed_u, &kwind.closed_u );
  w = xge_NewSwitch ( 1, w, sw13CLOSED_V, 109, 16, 0, 102,
                      txtClosed_v, &kwind.closed_v );
  w = xge_NewSwitch ( 1, w, sw13DOMAIN_NET, 109, 16, 0, 122,
                      txtDomain_net, &display_domain_net );
  w = xge_NewSwitch ( 1, w, sw13MARK_UNMARK, 109, 16, 0, 142,
                      txtMark_unmark, &kwind.selecting_mode );
  w = xge_NewButton ( 1, w, btn13EQUIDIST_U, 82, 19, 0, 162,
                      txtEquidist_u );
  w = xge_NewButton ( 1, w, btn13EQUIDIST_V, 82, 19, 0, 183,
                      txtEquidist_v );
  w = xge_NewButton ( 1, w, btn13FLIP, 82, 19, 0, 204, txtFlip );
  w = xge_NewSwitch ( 1, w, sw13COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &kwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 1, w, sw13PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &kwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 1, w, sw13MOVE_MANY_KNOTS, 109, 16, 0, xge_HEIGHT-76,
                      txtMove_many_knots, &kwind.moving_many );
  xge_SetWidgetPositioning ( w, 2, 0, -60 );
  menu4 = xge_NewMenu ( 1, menu3, MENU4, 110, xge_HEIGHT-36, 0, 20, w );

  w = xge_NewButton ( 1, NULL, btn16GENERAL, 76, 19, 2, 22, txtGeneral );
  w = xge_NewButton ( 1, w, btn16SPHERICAL, 76, 19, 2, 42, txtSpherical );
  w = xge_NewButton ( 1, w, btn16SWEPT, 76, 19, 2, 62, txtSwept );
  w = xge_NewButton ( 1, w, btn16BLENDING, 76, 19, 2, 82, txtBlending );
  popup10 = xge_NewFMenu ( 1, NULL, POPUP10, 80, 83, 0, 20, w );
  popup10->msgproc = xge_PopupMenuMsg;


  /* the first alternative for the second window - spherical product    */
  /* here we have two curve windows and two knot windows, attached      */
  /* alternatively to the widgets; one pair is for the equator, and     */
  /* the other os for the meridian. Initially the mappings are the same */
  eq_cwind.er = xge_New2Dwind ( 1, NULL, 1, xge_WIDTH-110, xge_HEIGHT-78,
                     110, 20, &eq_cwind, RysujSPROkno );
  eq_ckwind.er = xge_NewKnotWind ( 1, eq_cwind.er, 1, xge_WIDTH-110, 56,
                     110, xge_HEIGHT-56, &eq_ckwind, MAX_KNOTS, equator_knots );
  memcpy ( &mer_cwind, &eq_cwind, sizeof(xge_2Dwind) );
  memcpy ( &mer_ckwind, &eq_ckwind, sizeof(xge_KnotWind) );
  mer_ckwind.knots = &meridian_knots[0];

    /* the top menu */
  w = xge_NewButton ( 1, NULL, btn14bTYPE, 76, 19, 0, 0, txtSurface_Type );
  w = xge_NewButton ( 1, w, btn14bEDIT, 76, 19,  78, 0, txtEdit );
  w = xge_NewButton ( 1, w, btn14bVIEW, 76, 19, 156, 0, txtView );
  w = xge_NewButton ( 1, w, btn14bARCS, 76, 19, 234, 0, txtArcs );
  menu5 = xge_NewMenu ( 1, eq_ckwind.er, 9, xge_WIDTH, 20, 0, 0, w );
    /* the side menu - Edit the curves for spherical product surfaces */
  w = xge_NewTextWidget ( 1, NULL, txt15aSPHERICAL, 109, 16, 0, 22,
                          txtSpherical_product );
  w = xge_NewSwitch ( 1, w, sw15aBIND, 109, 16, 0, 40, txtBind, &bind_spr );
  w = xge_NewSwitch ( 1, w, sw15aEQUATOR, 109, 16, 0, 60, txtEquator, &equator );
  w = xge_NewSwitch ( 1, w, sw15aMERIDIAN, 109, 16, 0, 80, txtMeridian, &meridian );
  w = xge_NewButton ( 1, w, btn15aRESET, 76, 19, 0, 100, txtReset );
  w = xge_NewIntWidget ( 1, w, intw15aDEGREE, 76, 19, 0, 122,
                         0, MAX_DEGREE+1, &eqdeg, txtDegree, &eq_ckwind.degree );
  w = xge_NewSwitch ( 1, w, sw15aCLOSED, 109, 16, 0, 142,
                      txtClosed, &eqmer_closed );
  w = xge_NewSwitch ( 1, w, sw15aNURBS, 109, 16, 0, 162,
                      txtNURBS, &eqmer_nurbs );
  w = xge_NewSwitch ( 1, w, sw15aMARK_UNMARK, 109, 16, 0, 182,
                      txtMark_unmark, &eq_cwind.selecting_mode );
  w = xge_NewSwitch ( 1, w, sw15aMOVE, 109, 16, 0, 202,
                      txtMove, &eq_cwind.moving_tool );
  w = xge_NewSwitch ( 1, w, sw15aSCALE, 109, 16, 0, 222,
                      txtScale, &eq_cwind.scaling_tool );
  w = xge_NewSwitch ( 1, w, sw15aROTATE, 109, 16, 0, 242,
                      txtRotate, &eq_cwind.rotating_tool );
  w = xge_NewSwitch ( 1, w, sw15aSHEAR, 109, 16, 0, 262,
                      txtShear, &eq_cwind.shear_tool );
  w = xge_NewSwitch ( 1, w, sw15aCOORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &eq_cwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 1, w, sw15aPAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &eq_cwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 1, w, sw15aMOVE_MANY_KNOTS, 109, 16, 0, xge_HEIGHT-76,
                      txtMove_many_knots, &eq_ckwind.moving_many );
  xge_SetWidgetPositioning ( w, 2, 0, -60 );
  menu60list = w;
    /* the side menu - View the curves for spherical product surfaces */
  w = xge_NewTextWidget ( 1, NULL, txt15aSPHERICAL, 109, 16, 0, 22,
                          txtSpherical_product );
  w = xge_NewSwitch ( 1, w, sw15bCONTROL_POLYGON, 109, 16, 0, 40,
                      txtControl_polygon, &display_eqmer_control_polygon );
  w = xge_NewSwitch ( 1, w, sw15bCURVE, 109, 16, 0, 60,
                      txtCurve, &display_eqmer_curve );
  w = xge_NewSwitch ( 1, w, sw15bBEZIER_POLYGONS, 109, 16, 0, 80,
                      txtBezier_polygons, &display_eqmer_Bezier_polygons );
  w = xge_NewSwitch ( 1, w, sw15bTICKS, 109, 16, 0, 80,
                      txtTicks, &display_eqmer_ticks );
  w = xge_NewSwitch ( 1, w, sw15aCOORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &eq_cwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 1, w, sw15aPAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &eq_cwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 1, w, sw15aMOVE_MANY_KNOTS, 109, 16, 0, xge_HEIGHT-76,
                      txtMove_many_knots, &eq_ckwind.moving_many );
  xge_SetWidgetPositioning ( w, 2, 0, -60 );
  menu61list = w;
    /* the side menu - Arcs */
  w = xge_NewTextWidget ( 1, NULL, txt15aSPHERICAL, 109, 16, 0, 22,
                          txtSpherical_product );
  w = xge_NewButton ( 1, w, btn15cQUARTER_CIRCLE, 76, 19, 0, 40,
                      txtQuarter );
  w = xge_NewButton ( 1, w, btn15cHALF_CIRCLE, 76, 19, 0, 60,
                      txtHalf );
  w = xge_NewButton ( 1, w, btn15cFULL_CIRCLE, 76, 19, 0, 80,
                      txtFull_circle );
  w = xge_NewDiald ( 1, w, dial15cARC_ANGLE, 109, 39, 0, 100,
                     txtArc_angle, &arc_angle );
  w = xge_NewSwitch ( 1, w, sw15aCOORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &eq_cwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 1, w, sw15aPAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &eq_cwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 1, w, sw15aMOVE_MANY_KNOTS, 109, 16, 0, xge_HEIGHT-76,
                      txtMove_many_knots, &eq_ckwind.moving_many );
  xge_SetWidgetPositioning ( w, 2, 0, -60 );
  menu62list = w;

  menu6 = xge_NewMenu ( 1, menu5, MENU6, 110, xge_HEIGHT-36, 0, 20, menu60list );

  /* the second alternative for the second window --- swept surfaces */
    /* the top menu */
  w = xge_NewButton ( 1, NULL, btn14aTYPE, 76, 19, 0, 0, txtSurface_Type );
  menu9 = xge_NewMenu ( 1, NULL, MENU9, xge_WIDTH, 20, 0, 0, w );
    /* the side menu */
  w = xge_NewTextWidget ( 1, NULL, txt110SWEPT, 109, 16, 0, 22,
                          txtSwept_Surface );
  w = xge_NewSwitch ( 1, w, sw13COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &kwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 1, w, sw13PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &kwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  menu10 = xge_NewMenu ( 1, menu9, MENU10, 110, xge_HEIGHT-36, 0, 20, w );

  /* the third alternative for the second window --- blending surfaces */
    /* the top menu */
  w = xge_NewButton ( 1, NULL, btn14aTYPE, 76, 19, 0, 0, txtSurface_Type );
  w = xge_NewButton ( 1, w, btn111TRANSFORM, 76, 19, 78, 0, txtTransform );
  w = xge_NewButton ( 1, w, btn111OPTIONS, 76, 19, 156, 0, txtOptions );
  w = xge_NewButton ( 1, w, btn111INFO, 76, 19, xge_WIDTH-76, 0, txtInfo );
  xge_SetWidgetPositioning ( w, 1, -76, 0 );
  menu11 = xge_NewMenu ( 1, NULL, MENU11, xge_WIDTH, 20, 0, 0, w );
    /* the side menu */
  w = xge_NewTextWidget ( 1, NULL, txt112BLENDING, 109, 16, 0, 22,
                          txtBlending_Surface );
  w = xge_NewSwitch ( 1, w, sw112BLENDING_G1, 36, 16, 0, 40, txtG1, &sw_blending_g1 );
  w = xge_NewSwitch ( 1, w, sw112BLENDING_G2, 36, 16, 36, 40, txtG2, &sw_blending_g2 );
  w = xge_NewButton ( 1, w, btn112INIT_BLENDING, 76, 19, 0, 60, txtInit );
  w = xge_NewButton ( 1, w, btn112BLENDING_REFINE, 76, 19, 0, 81, txtRefine );
  w = xge_NewSwitch ( 1, w, sw112CLAMPED_BLENDING, 109, 16, 0, 102,
                      txtClamped, &sw_clamped_blending );
  w = xge_NewSwitch ( 1, w, sw112CLOSED_BLENDING, 76, 16, 0, 122,
                      txtClosed, &kwind.closed_u );
  w = xge_NewSwitch ( 1, w, sw112CONSTRAINTS, 109, 16, 0, 142,
                      txtConstraints, &sw_blending_constraints );
  w = bl_trihsw =
      xge_NewSwitch ( 1, w, sw112TRIHARMONIC_BLENDING, 109, 16, 0, 162,
                      txtTriharmonic, &sw_triharmonic_blending );
  w = xge_NewSlidebard ( 1, w, sl112NONLIN_BLENDING_C, 109, 10, 0, 182,
                         &blending_factor );
  w = bl_optbtn =
      xge_NewButton ( 1, w, btn112BLENDING_LMT_ITER, 76, 19, 0, 196, txtOptimize );
  w = xge_NewSwitch ( 1, w, sw112BLENDING_OPT_ENTIRE, 16, 16, 29, 239,
                      NULL, &sw_blending_opt_entire );
  w = xge_NewIntWidget ( 1, w, intw112BLENDING_UMIN, 23, 19, 0, 238,
                         3, MAX_KNOTS, &blending_opt_range[0], NULL,
                         &blending_opt_part[0] );
  w = xge_NewIntWidget ( 1, w, intw112BLENDING_UMAX, 23, 19, 50, 238,
                         3, MAX_KNOTS, &blending_opt_range[1], NULL,
                         &blending_opt_part[1] );
  w = xge_NewIntWidget ( 1, w, intw112BLENDING_VMIN, 23, 19, 25, 258,
                         3, MAX_KNOTS, &blending_opt_range[2], NULL,
                         &blending_opt_part[2] );
  w = xge_NewIntWidget ( 1, w, intw112BLENDING_VMAX, 23, 19, 25, 218,
                         3, MAX_KNOTS, &blending_opt_range[3], NULL,
                         &blending_opt_part[3] );
  w = xge_NewSwitch ( 1, w, sw13COORDINATES, 109, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &kwind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -20 );
  w = xge_NewSwitch ( 1, w, sw13PAN_ZOOM, 109, 16, 0, xge_HEIGHT-56,
                      txtPan_zoom, &kwind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -40 );
  w = xge_NewSwitch ( 1, w, sw13MARK_UNMARK, 109, 16, 0, 76,
                      txtMark_unmark, &kwind.selecting_mode );
  xge_SetWidgetPositioning ( w, 2, 0, -60 );
  menu12 = xge_NewMenu ( 1, menu11, MENU12, 110, xge_HEIGHT-36, 0, 20, w );

    /* the bottom menu */
  status1sw = xge_NewSwitch ( 1, NULL, sw16STATUS, 109, 16, 0, xge_HEIGHT-16,
                              txtNULL, &status_1_sw );
  xge_SetWidgetPositioning ( status1sw, 0, 0, 0 );
  status1 = xge_NewTextWidget ( 1, status1sw, STATUSLINE1, xge_MAX_WIDTH-110, 14,
                                110, xge_HEIGHT-14, statustext1 );
  xge_SetWidgetPositioning ( status1, 0, 110, 1 );
  menu7 = xge_NewMenu ( 1, menu4, MENU7, xge_WIDTH, 16, 0, xge_HEIGHT-16, status1 );
    /* a copy for the alternative window contents */
  w = xge_NewSwitch ( 1, NULL, sw16STATUS, 109, 16, 0, xge_HEIGHT-16,
                      txtNULL, &status_1_sw );
  xge_SetWidgetPositioning ( w, 0, 0, 0 );
  menu8 = xge_NewMenu ( 1, menu6, MENU8, xge_WIDTH, 16, 0, xge_HEIGHT-16, w );
    /* a copy for the second alternative window contents */
  w = xge_NewSwitch ( 1, NULL, sw16STATUS, 109, 16, 0, xge_HEIGHT-16,
                      txtNULL, &status_1_sw );
  xge_SetWidgetPositioning ( w, 0, 0, 0 );
  menu13 = xge_NewMenu ( 1, menu10, MENU13, xge_WIDTH, 16, 0, xge_HEIGHT-16, w );
    /* a copy for the third alternative window contents */
  w = xge_NewSwitch ( 1, NULL, sw16STATUS, 109, 16, 0, xge_HEIGHT-16,
                      txtNULL, &status_1_sw );
  xge_SetWidgetPositioning ( w, 0, 0, 0 );
  menu14 = xge_NewMenu ( 1, menu12, MENU14, xge_WIDTH, 16, 0, xge_HEIGHT-16, w );

    /* setup the Transform popup menu */  
  w = xge_NewStringEd ( 1, NULL, txtedP11A0, 74, 19, 20+164, 40+10,
                        12, blending_trans_str[0], &blending_trans_editor[0] );
  w = xge_NewStringEd ( 1, w, txtedP11A1, 74, 19, 20+242, 40+10,
                        12, blending_trans_str[1], &blending_trans_editor[1] );
  w = xge_NewStringEd ( 1, w, txtedP11A2, 74, 19, 20+320, 40+10,
                        12, blending_trans_str[2], &blending_trans_editor[2] );
  w = xge_NewStringEd ( 1, w, txtedP11A3, 74, 19, 20+164, 40+31,
                        12, blending_trans_str[3], &blending_trans_editor[3] );
  w = xge_NewStringEd ( 1, w, txtedP11A4, 74, 19, 20+242, 40+31,
                        12, blending_trans_str[4], &blending_trans_editor[4] );
  w = xge_NewStringEd ( 1, w, txtedP11A5, 74, 19, 20+320, 40+31,
                        12, blending_trans_str[5], &blending_trans_editor[5] );
  w = xge_NewStringEd ( 1, w, txtedP11A6, 74, 19, 20+164, 40+52,
                        12, blending_trans_str[6], &blending_trans_editor[6] );
  w = xge_NewStringEd ( 1, w, txtedP11A7, 74, 19, 20+242, 40+52,
                        12, blending_trans_str[7], &blending_trans_editor[7] );
  w = xge_NewStringEd ( 1, w, txtedP11A8, 74, 19, 20+320, 40+52,
                        12, blending_trans_str[8], &blending_trans_editor[8] );
  w = xge_NewButton ( 1, w, btnP11IDENTITY, 76, 19, 20+6, 40+10, txtIdentity );
  w = xge_NewButton ( 1, w, btnP11OK, 76, 19, 20+162, 40+180-30, txtOK );
  popup11 = xge_NewFMenu ( 1, NULL, POPUP11, 400, 180, 20, 40, w );

    /* setup the Options popup menu */
  w = xge_NewIntWidget ( 1, NULL, intwP12BLENDING_LMT_NITER, 115, 19, 20+4, 40+4,
                         1, 99, &wdg_blending_lmt_iter, txtIterationLimit,
                         &blending_lmt_iter );
  w = xge_NewIntWidget ( 1, w, intwP12BLENDING_QUAD1, 115, 19, 20+4, 40+24,
                 3, 10, &wdg_blending_quad1, txtQuadratureKnots, &blending_quad1 );
  w = xge_NewIntWidget ( 1, w, intwP12BLENDING_QUAD2, 17, 19, 20+4+115+2, 40+24,
                         3, 10, &wdg_blending_quad2, txtNULL, &blending_quad2 );
  w = xge_NewSwitch ( 1, w, swP12SHOWSTEPS, 109, 16, 20+4, 40+45, txtShowSteps,
                      &sw_show_steps );
  w = xge_NewSwitch ( 1, w, swP12DUMPDATA, 109, 16, 20+4, 40+65, txtDumpData,
                      &sw_blending_opt_dump );
  w = xge_NewButton ( 1, w, btnP12OK, 76, 19, 20+33, 40+86, txtOK );
  popup12 = xge_NewFMenu ( 1, NULL, POPUP12, 142, 108, 20, 40, w );


  xge_SetWinEdRect ( menu7 );
  ResizeWinStatus ( win1 );
} /*SetupWindow1Widgets*/

void init_edwin ( void )
{
  int id;

  setvbuf ( stdout, NULL, _IONBF, 0 );
  if ( !pkv_InitScratchMem ( SCRATCH_MEM_SIZE ) ) {
    printf ( "Error: cannot allocate scratch memory stack\n" );
    exit ( 1 );
  }
  SetupWindow0Widgets ();
  SetupWindow1Widgets ();

  ResetObject ();
  SetKWindNKnots ();
  xge_T2KnotWindResetMapping ( &kwind );
  SelectEquator ();
  FindBoundingBox ( &swind.DefBBox );
  xge_3DwindInitProjections ( &swind, swind.DefBBox.x0, swind.DefBBox.x1,
      swind.DefBBox.y0, swind.DefBBox.y1, swind.DefBBox.z0, swind.DefBBox.z1 );
  for ( id = 0; id < 4; id++ )
    ProjectSurface ( id );

  RendInit ();
  xge_RedrawAll ();
  xge_SetWindow ( win0 );
  xge_DisplayInfoMessage ( InfoMsg, -1 );
} /*init_edwin*/

void destroy_edwin ( void )
{
  void kill ( pid_t pid, int sig );

  RendDestroy ();
  printf ( "Scratch memory used: %d out of %d bytes\n",
           (int)pkv_MaxScratchTaken(), SCRATCH_MEM_SIZE );
  pkv_DestroyScratchMem ();
  if ( xge_ChildIsActive () ) {
    signal ( SIGCHLD, SIG_DFL );    
    kill ( xge_child_pid, SIGKILL );
  }
} /*destroy_edwin*/

