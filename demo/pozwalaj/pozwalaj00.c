
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "bsmesh.h"
#include "g2blendingd.h"
#include "egholed.h"
#include "mengerc.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

#include "widgets.h"
#include "editor.h"
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"
#include "pozwalaj.h"

char **prog_argv = NULL;
int  ncpu = 1;  /* number of CPUs - changed by program initialisation */

int  win0, win1;

/* window 0 widgets and stuff */
xge_widget    *top00menu = NULL, *side00menu = NULL, *bottom00menu = NULL,
              *geom00menu = NULL, *geom00win = NULL,
              *popup00 = NULL, *popup01 = NULL, *popup02 = NULL,
              *popup03 = NULL;
boolean       side00menu_wide = false;

xge_listbox   dirlist1, filelist1, dirlist2, filelist2;
xge_string_ed filename_editor;
const char    file_filter[] = "*.bs";
const char    file_ext[] = ".bs";
char          initial_directory[MAX_PATH_LGT+1],
              current_directory[MAX_PATH_LGT+1],
              current_dir[MAX_PATH_SHRT+1];
char          filename[MAX_FILENAME_LGT+1];
boolean       showhiddenfiles = false;

xge_2Dwind    g00win2D;
xge_widget    *geom00win2D = NULL;
xge_3Dwind    g00win3D;
xge_widget    *geom00win3D = NULL;

boolean       markbits[5] = {true,false,false,false,false},
              swwin0mark = false, swwin0translate = false, swwin0scale = false,
              swwin0rotate = false, swwin0shear = false, swwin0panzoom = false,
              swwin0coordinates = false,
              sw_bsm_selectvertex, sw_bsm_selectedge;

xge_widget    *win0statl = NULL, *win0cmdl = NULL;
boolean       win0statusline = false, win0commandline = false;
xge_string_ed command0_editor;
char          status0[MAX_COMMAND_LGT+1], command0[MAX_COMMAND_LGT+1];

xge_widget     *side00widgets = NULL,
               *side00awidgets = NULL, *side00dwidgets = NULL,
               *side00escroll = NULL, *side00econtents = NULL, *side00ewidgets = NULL;
xge_scroll_widget side00esw;

xge_string_ed  side00eparam_ed[8];
char           side00eparam_str[8][MAX_PARAM_LGT+1] =
               {"0.0", "0.0", "0.0", "1.0", "0.0", "0.0", "1.0", "0.0'0\""};
double         side00eparam[8] = {0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0};

  /* rendering menu stuff */
xge_widget     *side00bwidgets = NULL, *side00cwidgets = NULL, *side00fwidgets = NULL;
xge_widget     *renderbtn0 = NULL, *renderbtn1 = NULL, *renderbtn2 = NULL;
xge_int_widget side00bnpthreads;
int            rendering_npthreads = 1;

  /* camera parameters */
double         side00fpsi = 0.0, side00ftheta = 0.0, side00fphi = 0.0;
xge_string_ed  side00fparam_ed[7];
char           side00fparam_str[7][MAX_PARAM_LGT+1] =
               { "", "", "", "", "", "", "" };
  /* file saving options */
boolean        sw_save_all = true, sw_save_active = false, sw_save_current = false,
               sw_save_camera = false, sw_save_append = false;

/* window 1 widgets and stuff */
char           whichside10menu = SIDE10MENU_EDIT;
boolean        side10menu_wide = false;

xge_widget      *top10menu = NULL, *side10menu = NULL, *bottom10menu = NULL,
                *geom10menu = NULL, *geom10win = NULL,
                *popup10 = NULL, *popup11 = NULL, *popup12 = NULL, *popup13 = NULL,
                *popup14 = NULL, *popup15 = NULL;

xge_widget      *side10wdg_none = NULL,
                *side10wdg_bezc_edit = NULL, *side10wdg_bsc_edit = NULL,
                *side10wdg_bezp_edit = NULL, *side10wdg_bsp_edit = NULL,
                *side10wdg_bsm_editscroll = NULL, *side10wdg_bsm_editcontents = NULL,
                *side10wdg_bsm_edit = NULL,
                *side10wdg_bsh_edit = NULL,
                *geom10wdg_bezc = NULL, *geom10wdg_bsc = NULL, *geom10wdg_bezp = NULL,
                *geom10wdg_bsp = NULL, *geom10wdg_bsm = NULL, *geom10wdg_bsh = NULL;
xge_scroll_widget side10_bsm_editsw, side10_bsm_optsw;
xge_widget      *side10wdg_bezc_view = NULL, *side10wdg_bsc_view = NULL,
                *side10wdg_bezp_view = NULL, *side10wdg_bsp_view = NULL,
                *side10wdg_bsm_view = NULL, *side10wdg_bsh_view = NULL;
xge_widget      *side10wdg_bezc_data = NULL, *side10wdg_bsc_data = NULL,
                *side10wdg_bezp_data = NULL, *side10wdg_bsp_data = NULL,
                *side10wdg_bsm_data = NULL, *side10wdg_bsh_data = NULL;
xge_widget      *side10wdg_bezc_opt = NULL,
                *side10wdg_bsc_opt = NULL,
                *side10wdg_bsc_opt1 = NULL, *side10wdg_bsc_opt2 = NULL,
                *side10wdg_bezp_opt = NULL,
                *side10wdg_bsp_opt = NULL,
                *side10wdg_bsp_opt_general = NULL,
                *side10wdg_bsp_opt_spherical = NULL,
                *side10wdg_bsp_opt_blendingG1 = NULL,
                *side10wdg_bsp_opt_blendingG2 = NULL,
                *side10wdg_bsm_optscroll = NULL,
                *side10wdg_bsm_optcontents = NULL,
                *side10wdg_bsm_opt = NULL, *side10wdg_bsm_opt1 = NULL,
                *side10wdg_bsm_opt2 = NULL, *side10wdg_bsm_opt3 = NULL,
                *side10wdg_bsh_opt = NULL;

xge_2Dwind      g10win2D, g10win2Deqmer;
xge_widget      *geom10win2D = NULL, *geom10win2Deqmer = NULL;
xge_KnotWind    g10knotwin, g10knotwineqmer;
xge_widget      *geom10knotwin = NULL, *geom10knotwineqmer = NULL;
xge_T2KnotWind  g10t2knotwin;
xge_widget      *geom10t2knotwin = NULL;

xge_widget      *popup11wdg_bezc = NULL, *popup11wdg_bezp = NULL,
                *popup11wdg_bsc = NULL, *popup11wdg_bsp = NULL,
                *popup11wdg_bsm = NULL, *popup11wdg_bsh = NULL;

xge_widget      *win1statl = NULL, *win1cmdl = NULL;
boolean         win1statusline = false, win1commandline = false;
xge_string_ed   command1_editor;
char            status1[MAX_COMMAND_LGT+1], command1[MAX_COMMAND_LGT+1];

xge_string_ed   objname_editor1, objname_editor2;
char            objectname[MAX_NAME_LENGTH+1];

xge_listbox     objectlist;

xge_int_widget  intwdensityu, intwdensityv;
int             density_u, density_v;

  /* object type specific stuff */
int             degree, degreeu, degreev;
boolean         sw_view_curve = true, sw_view_cpoly = true,
                sw_view_surf = true, sw_view_cnet = true,
                sw_view_curvature = false, sw_view_torsion = false;
int             curv_graph_dens = 10;
double          curvature_scale = 0.0, torsion_scale = 0.0;
boolean         sw_view_hole_filling = true,
                sw_hfill_g1 = true, sw_hfill_g2 = false, sw_hfill_g1q2 = false;
boolean         sw_hfill_coons = true, sw_hfill_bezier = false;
double          sl_g1q2_param = 0.5;
double          sl_pipe_diameter = 0.5;

    /* Bezier curves */
xge_string_ed   bezc_name_ed;
xge_int_widget  intw_bezcdeg;

    /* Bezier patches */
xge_string_ed   bezp_name_ed;
xge_int_widget  intw_bezpdegu, intw_bezpdegv;
xge_widget      *bezp_dens_u, *bezp_dens_v;

    /* B-spline curves */
xge_string_ed   bsc_name_ed;
xge_int_widget  intw_bscdeg, bsc_graphdens;
boolean         bsc_sw_closed = false;
boolean         bsc_sw_view_bpoly = false;
boolean         bsc_sw_mengerc = false,
                bsc_sw_mc_logit = false;
double          bsc_sl_mcexp, bsc_sl_mcppar[5];
xge_int_widget  bsc_mc_qkn, bsc_mc_popt, bsc_mc_maxit, bsc_mc_nthr,
                bscnpthreads;
int             bsc_mc_qknots = 4, bsc_mc_ppopt = 3, bsc_mc_maxiter = 50,
                bsc_npthreads;

    /* B-spline patches */
xge_string_ed   bsp_name_ed;
xge_int_widget  intw_bspdegu, intw_bspdegv;
xge_widget      *bsp_dens_u = NULL, *bsp_dens_v = NULL;
boolean         bsp_sw_closed_u = false, bsp_sw_closed_v = false,
                sw_bsp_clamped = false;
boolean         sw_bsp_dom_coord = false, sw_bsp_dom_panzoom = false;
boolean         bsp_bl_nharmonic = false, bsp_blending_entire = true;
xge_widget      *sw_bsp_bl_nharmonic = NULL;
xge_int_widget  bsp_blG1_minu, bsp_blG1_maxu, bsp_blG1_minv, bsp_blG1_maxv,
                bsp_blG2_minu, bsp_blG2_maxu, bsp_blG2_minv, bsp_blG2_maxv;
int             bsp_bl_opt_range[4];
double          sl_bsp_bl_param = 0.5;
int             bsp_bl_nkn1 = 6, bsp_bl_nkn2 = 8, bsp_bl_maxit = 20;
xge_int_widget  bsp_bl_qknots1, bsp_bl_qknots2, bsp_bl_maxiter;
xge_widget      *bsp_bl_optimizeG1 = NULL, *bsp_bl_optimizeG2 = NULL;
xge_widget      *bsp_type_button = NULL;

boolean         sw_bsp_sproduct_equator     = true,
                sw_bsp_sproduct_meridian    = false,
                sw_bsp_sproduct_eqrational  = false,
                sw_bsp_sproduct_eqclosed    = false,
                sw_bsp_sproduct_equniform   = false,
                sw_bsp_sproduct_merrational = false,
                sw_bsp_sproduct_merclosed   = false,
                sw_bsp_sproduct_meruniform  = false;
xge_string_ed   bsp_sproduct_eqname_ed, bsp_sproduct_mername_ed;
xge_int_widget  intw_bsp_sproduct_eqdeg, intw_bsp_sproduct_merdeg;
int             bsp_eqdegree = 1, bsp_merdegree = 1;
GO_BSplineCurve *bsp_sproduct_eqmer;

    /* B-spline meshes */
xge_string_ed   bsm_name_ed;
xge_int_widget  intw_bsmdeg;
xge_int_widget  intw_bsm_vert0, intw_bsm_vert1,
                intw_bsm_edge0, intw_bsm_edge1,
                intw_bsm_fac0, intw_bsm_fac1;
int             bsm_vertex_num0 = -1, bsm_vertex_num1 = -1,
                bsm_edge_num0 = -1, bsm_edge_num1 = -1,
                bsm_facet_num0 = -1, bsm_facet_num1 = -1;
xge_int_widget  intw_bsmdeg1, intw_bsmdeg2;
int             bsm_degree1 = 3, bsm_degree2 = 3;
xge_widget      *bsm_density;

boolean         bsm_sw_subdivision = false;
boolean         bsm_sw_blending = false;
boolean         bsm_sw_view_special = false;

boolean         bsm_sw_constr = false, bsm_sw_log_it = false,
                bsm_sw_shape_only = false,
                bsm_sw_alt_multilevel = false;
double          sl_bsm_bl_param = 0.5;
int             bsm_bl_nkn1 = 6, bsm_bl_nkn2 = 8, bsm_bl_maxit = 20,
                bsm_bl_nlev = 0, bsm_bl_nbl = 1;
xge_int_widget  bsm_bl_qknots1, bsm_bl_qknots2, bsm_bl_maxiter,
                bsm_bl_nlevels, bsm_bl_nblocks;
xge_widget      *bsm_bl_optimize, *bsc_mc_optimize;
boolean         bsm_sw_use_coarse = false;
xge_string_ed   bsm_coarse_editor;
xge_int_widget  bsmnpthreads; 
int             bsm_npthreads = 1;

boolean         bsm_sw_data_add = false, bsm_sw_data_replace = true;

char            logfilename[] = "pozwalaj.log.bs";

    /* B-spline holes */

    /* colour specimen widget */
xge_widget     *colour_box = NULL;
double         colour_rgb[3] = {1.0,1.0,0.0};

    /* ray traced picture data */
boolean rendered_picture = false;
    /* lights & special stuff editing mode */
boolean picture_lights = true, editing_shapefunc = false, editing_lights = false;
boolean editing_pretrans = false, editing_camera = false;


    /* popup menu 14 - pre-transformation */
xge_string_ed pretrans_ed[12];
char          pretrans_str[12][MAX_PARAM_LGT+1];

