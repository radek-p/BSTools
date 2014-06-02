
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

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "xgedit.h"

#include "render.h"
#include "spl3d.h"
#include "ed3ds.h"
#include "pomnijipc.h"

/* ///////////////////////////////////////////////////////////////////////// */
int  prog_argc;
char **prog_argv;

int win0, win1;  /* window identifiers */

xge_widget *menu0, *menu1;
xge_3Dwind swind;
xge_widget *menu00list, *menu01list, *menu02list, *menu03list;
xge_widget *popup00, *popup01, *popup02, *popup03, *popup04;
xge_int_widget ddeg_u, ddeg_v, ddens_u, ddens_v;

xge_widget *domwind, *menu2, *menu3;
xge_T2KnotWind kwind;

xge_2Dwind     eq_cwind, mer_cwind;
xge_KnotWind   eq_ckwind, mer_ckwind;
xge_int_widget eqdeg;

xge_widget *menu4, *menu5, *menu6, *menu7, *menu8, *menu9,
           *menu10, *menu11, *menu12, *menu13, *menu14;

xge_widget *menu60list, *menu61list, *menu62list;
xge_widget *popup10, *popup11, *popup12;

xge_listbox filelist, dirlist;
xge_string_ed filename_editor;

/* ///////////////////////////////////////////////////////////////////////// */
const char file_filter[] = "*.bs";
const char file_ext[] = ".bs";

char current_directory[MAX_PATH_LGT+1] = "qq";
char filename[MAX_FILENAME_LGT+1] = "";

/* ///////////////////////////////////////////////////////////////////////// */
xge_int_widget wdg_blending_lmt_iter;
xge_widget *bl_optbtn;
int blending_quad1 = 6, blending_quad2 = 8;
xge_int_widget wdg_blending_quad1, wdg_blending_quad2;

xge_widget *bl_trihsw;
xge_int_widget blending_opt_range[4];
xge_string_ed blending_trans_editor[9];
char blending_trans_str[9][13] =
  {"1.0","0.0","0.0","0.0","1.0","0.0","0.0","0.0","1.0"};
boolean sw_show_steps = true;

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *status0, *status0sw;
boolean status_0_sw = false;
char statustext0[STATUS_LINE_LENGTH+1]; 

xge_widget *status1, *status1sw;
boolean status_1_sw = false;
char statustext1[STATUS_LINE_LENGTH+1];

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *renderbtn0, *renderbtn1;

boolean edit_light_sw[4] = {false,false,false,false};
boolean edit_reflection_frame = false;
boolean edit_highlight_frame = false;
boolean edit_sections_frame = false;

boolean swShadows = false, swAntialias = false;
boolean swind_picture = false;

double dfs_factor = 0.5;
double cf_range[2] = {0.0, 1.0};

/* light directions and intensities - the last is the ambient intensity */
vector3d light_dir[R_NLIGHTS] = {{0.0,0.0,1.0},
                                 {1.0,0.0,0.0},
                                 {0.0,1.0,0.0},
                                 {0.0,0.0,-1.0}};
double   light_int[R_NLIGHTS+1] = {0.9,0.0,0.0,0.0,0.15};

/* ///////////////////////////////////////////////////////////////////////// */
#define PATH_MAX 512
boolean child_ready = false;
static char initial_dir[PATH_MAX];

boolean MakeTheChildProcess ( void )
{
  child_ready = false;
  if ( chdir ( initial_dir ) )
    return false;
  return xge_MakeTheChild ( prog_argv[0], "_proc", POMNIJ_IPC_MAGIC );
} /*MakeTheChildProcess*/

int main ( int argc, char *argv[] )
{
  prog_argc = argc;
  prog_argv = argv;
  xge_Init ( argc, argv, CallBack, NULL );
  if ( !getcwd ( initial_dir, PATH_MAX ) ) {
    printf ( "Failed to get initial directory\n" );
    exit ( 1 );
  }
  if ( !MakeTheChildProcess () ) {
    xge_Cleanup ();
    exit ( 1 );
  }
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();
  xge_Cleanup ();
  exit ( 0 );
} /*main*/

