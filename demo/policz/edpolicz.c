
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

int win0, win1;
xge_widget *menu0, *menu1, *menu2;
xge_widget *menu01alist, *menu01blist, *menu01clist,
           *menu01dlist, *menu01elist, *menu01flist;
xge_widget *popup00, *popup01, *popup02, *popup03, *popup04;
xge_widget *status0, *status0sw, *renderbtn0, *renderbtn1;
xge_int_widget hole_k_intw;

/* constraint switches - to be reconfigured whenever necessary */
xge_widget *menu01cscrolled;
xge_scroll_widget scroll_constr_sw;
xge_widget *constr_sw[NUM_CONSTR_SWITCHES];
boolean constraint[NUM_CONSTR_SWITCHES];   
xge_widget *constr_row_text[MAX_CONSTR_SWITCH_ROWS];
boolean constraints1 = true;      /* first surface */
boolean constraints2 = false;     /* second surface */
boolean constraints1st = false;   /* first type */
boolean constraints2nd = false;   /* second type */
boolean constraintsZero = false;

xge_widget *menu3, *menu4, *menu5;
xge_widget *menu14alist, *menu14a1list, *menu14a2list,
           *menu14blist, *menu14clist, *menu14dlist;
xge_widget *knotwind;
xge_widget *status1, *status1sw;

xge_listbox   filelist, dirlist;
xge_string_ed filename_editor;
char current_directory[MAX_PATH_LGT+1] = "";
char filename[MAX_FILENAME_LGT+4] = "";
const char file_ext[] = ".bs";
const char file_filter[] = "*.bs";

boolean status_0_sw = false;
boolean status_1_sw = false;

char statustext0[STATUS_LINE_LENGTH+1] = "";
char statustext1[STATUS_LINE_LENGTH+1] = "";

boolean sw_opt_1 = true;
boolean sw_opt_2 = false;
xge_int_widget opt_order_1, opt_nk_1, opt_m1_1, opt_m2_1;
xge_int_widget opt_order_2, opt_nk_2, opt_m1_2, opt_m2_2;

boolean view_constraints_frame = false;
boolean view_surf_cp = true;
boolean view_surf_spatches = true;
boolean view_surf_1 = false;
boolean view_surf_2 = false;
boolean view_surf_numbers = false;

boolean view_dom_cp = true;
boolean view_dom_spatches = true;
boolean view_dom_patches1 = true;
boolean view_dom_patches2 = false;
boolean view_dom_numbers = false;

boolean edit_light_sw[4] = {false,false,false,false};
boolean edit_reflection_frame = false;
boolean edit_highlight_frame = false;
boolean edit_sections_frame = false;

boolean swShadows = false, swAntialias = false;

/* light directions and intensities - the last is the ambient intensity */
vector3d light_dir[R_NLIGHTS] = {{0.0,0.0,1.0},
                                 {1.0,0.0,0.0},
                                 {0.0,1.0,0.0},
                                 {0.0,0.0,-1.0}};
double   light_int[R_NLIGHTS+1] = {1.0,0.0,0.0,0.0,0.15};

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

