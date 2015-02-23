
/* ///////////////////////////////////////////////////  */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* ///////////////////////////////////////////////////  */

#ifndef POZWALAJ0_H
#include "pozwalaj0.h"
#endif

#define MAX_COMMAND_LGT   xge_MAX_STRING_LENGTH
#define MAX_PARAM_LGT      20

#define TOPMENUHEIGHT      20
#define BOTTOMMENUHEIGHT0  16
#define BOTTOMMENUHEIGHT1  21
#define SIDEMENUWIDTH0    110
#define SIDEMENUWIDTH1    122

    /* identifiers of side menus in the domain (win1) window */
#define SIDE10MENU_EDIT     1
#define SIDE10MENU_VIEW     2
#define SIDE10MENU_DATA     3
#define SIDE10MENU_OPTIONS  4

    /* commands to be executed when the processor is idle */
#define IDLE_COMMAND_RENDER_START           1
#define IDLE_COMMAND_RENDER_CONT            2
#define IDLE_COMMAND_BSP_BLENDING_OPT_INIT  3
#define IDLE_COMMAND_BSM_BLENDING_OPT_INIT  4
#define IDLE_COMMAND_BSC_MENGERC_OPT_INIT   5
#define IDLE_COMMAND_SEND_DATA_TO_CHILD     6
#define IDLE_COMMAND_BSMESH_OPT_SPECIALS    7

extern int           win0, win1;

/* window 0 widgets and stuff */
extern xge_widget     *top00menu, *side00menu, *bottom00menu,
                      *geom00menu, *geom00win,
                      *popup00, *popup01, *popup02, *popup03;
extern boolean        side00menu_wide;

extern xge_listbox    dirlist1, filelist1, dirlist2, filelist2;
extern xge_string_ed  filename_editor;
extern const char     file_filter[], file_ext[];
extern char           initial_directory[], current_directory[], current_dir[];
extern char           filename[];
extern boolean        showhiddenfiles;

extern xge_2Dwind     g00win2D;
extern xge_widget     *geom00win2D;
extern xge_3Dwind     g00win3D;
extern xge_widget     *geom00win3D;

extern boolean        markbits[5],
                      swwin0mark, swwin0markedg, swwin0translate,
                      swwin0scale, swwin0rotate,
                      swwin0shear, swwin0panzoom, swwin0coordinates,
                      sw_bsm_selectvertex, sw_bsm_selectedge;

extern xge_widget     *win0statl, *win0cmdl;
extern boolean        win0statusline, win0commandline;
extern xge_string_ed  command0_editor;
extern char           status0[], command0[];

extern xge_widget     *side00widgets, *side00awidgets, *side00dwidgets,
                      *side00escroll, *side00econtents, *side00ewidgets;
extern xge_scroll_widget side00esw;

extern xge_string_ed  side00eparam_ed[8];
extern char           side00eparam_str[8][MAX_PARAM_LGT+1];
extern double         side00eparam[8];

  /* rendering menu stuff */
extern xge_widget     *side00bwidgets, *side00cwidgets, *side00fwidgets;
extern xge_widget     *renderbtn0, *renderbtn1, *renderbtn2;
extern xge_int_widget side00bnpthreads;
extern int            rendering_npthreads;
extern boolean        swGaussian_c, swMean_c, swLambert_c, swReflection_c,
                      swHighlight_c, swSections_c;
extern boolean        swGaussian_d, swMean_d, swLambert_d, swReflection_d,
                      swHighlight_d, swSections_d;
extern boolean        swAntialias, swShadows;
extern double         render_cfrange[2], render_dfsf;
extern pkRenderer     rend;
extern XImage         *rendimage;

  /* camera parameters */
extern double         side00fpsi, side00ftheta, side00fphi;
extern xge_string_ed  side00fparam_ed[7];
extern char           side00fparam_str[7][MAX_PARAM_LGT+1];

  /* file saving options */
extern boolean        sw_save_all, sw_save_active, sw_save_current,
                      sw_save_camera, sw_save_append;

/* window 1 widgets and stuff */
extern xge_widget     *top10menu, *side10menu, *bottom10menu, *geom10menu,
                      *geom10win,
                      *popup10, *popup11, *popup12, *popup13, *popup14, *popup15;

extern boolean        side10menu_wide;

extern char           whichside10menu;
extern xge_widget     *side10wdg_none,
                      *side10wdg_bezc_edit, *side10wdg_bsc_edit, *side10wdg_bezp_edit,
                      *side10wdg_bsp_edit,
                      *side10wdg_bsm_editscroll, *side10wdg_bsm_editcontents,
                      *side10wdg_bsm_edit,
                      *side10wdg_bsh_edit,
                      *geom10wdg_bezc, *geom10wdg_bsc, *geom10wdg_bezp,
                      *geom10wdg_bsp, *geom10wdg_bsm, *geom10wdg_bsh;
extern xge_scroll_widget side10_bsm_editsw, side10_bsm_optsw;
extern xge_widget     *side10wdg_bezc_view, *side10wdg_bsc_view, *side10wdg_bezp_view,
                      *side10wdg_bsp_view, *side10wdg_bsm_view, *side10wdg_bsh_view;
extern xge_widget     *side10wdg_bezc_data, *side10wdg_bsc_data, *side10wdg_bezp_data,
                      *side10wdg_bsp_data, *side10wdg_bsm_data, *side10wdg_bsh_data;
extern xge_widget     *side10wdg_bezc_opt,
                      *side10wdg_bsc_opt,
                      *side10wdg_bsc_opt1, *side10wdg_bsc_opt2,
                      *side10wdg_bezp_opt,
                      *side10wdg_bsp_opt, *side10wdg_bsp_opt_general,
                      *side10wdg_bsp_opt_spherical,
                      *side10wdg_bsp_opt_blendingG1,
                      *side10wdg_bsp_opt_blendingG2,
                      *side10wdg_bsm_optscroll, *side10wdg_bsm_optcontents,
                      *side10wdg_bsm_opt, *side10wdg_bsm_opt1, *side10wdg_bsm_opt2,
                      *side10wdg_bsm_opt3,
                      *side10wdg_bsh_opt;

extern xge_2Dwind     g10win2D, g10win2Deqmer;
extern xge_widget     *geom10win2D, *geom10win2Deqmer;
extern xge_KnotWind   g10knotwin, g10knotwineqmer;
extern xge_widget     *geom10knotwin, *geom10knotwineqmer;
extern xge_T2KnotWind g10t2knotwin;
extern xge_widget     *geom10t2knotwin;

extern xge_widget     *popup11wdg_bezc, *popup11wdg_bezp,
                      *popup11wdg_bsc, *popup11wdg_bsp,
                      *popup11wdg_bsm, *popup11wdg_bsh;

extern xge_widget     *win1statl, *win1cmdl;
extern boolean        win1statusline, win1commandline;
extern xge_string_ed  command1_editor;
extern char           status1[], command1[];

extern xge_string_ed  objname_editor1, objname_editor2;
extern char           objectname[];

extern xge_listbox    objectlist;

extern xge_int_widget intwdensityu, intwdensityv;
extern int            density_u, density_v;

  /* object type specific stuff */
extern int            degree, degreeu, degreev;
extern boolean        sw_view_curve, sw_view_cpoly, sw_view_surf, sw_view_cnet;
extern boolean        sw_view_curvature, sw_view_torsion;
extern double         curvature_scale, torsion_scale;
extern int            curv_graph_dens;
extern boolean        sw_view_hole_filling, sw_hfill_g1, sw_hfill_g2, sw_hfill_g1q2;
extern boolean        sw_hfill_coons, sw_hfill_bezier;
extern double         sl_g1q2_param;
extern double         sl_pipe_diameter;

    /* Bezier curves */
extern xge_string_ed  bezc_name_ed;
extern xge_int_widget intw_bezcdeg;

    /* Bezier patches */
extern xge_string_ed  bezp_name_ed;
extern xge_int_widget intw_bezpdegu, intw_bezpdegv;
extern xge_widget     *bezp_dens_u, *bezp_dens_v;

    /* B-spline curves */
extern xge_string_ed  bsc_name_ed;
extern xge_int_widget intw_bscdeg, bsc_graphdens;
extern boolean        bsc_sw_closed;
extern boolean        sw_bsc_move_many_knots, sw_bsc_dom_coord, sw_bsc_dom_panzoom;
extern boolean        bsc_sw_view_bpoly;
extern boolean        bsc_sw_mengerc, bsc_sw_mc_logit;
extern double         bsc_sl_mcexp, bsc_sl_mcppar[5];
extern xge_int_widget bsc_mc_qkn, bsc_mc_popt, bsc_mc_maxit, bsc_mc_nthr,
                      bscnpthreads;
extern int            bsc_mc_qknots, bsc_mc_ppopt, bsc_mc_maxiter,
                      bsc_npthreads;

    /* B-spline patches */
extern xge_string_ed  bsp_name_ed;
extern xge_int_widget intw_bspdegu, intw_bspdegv;
extern xge_widget     *bsp_dens_u, *bsp_dens_v;
extern boolean        bsp_sw_closed_u, bsp_sw_closed_v, sw_bsp_clamped;
extern boolean        sw_bsp_move_many_knots, sw_bsp_dom_coord, sw_bsp_dom_panzoom;
extern boolean        bsp_bl_nharmonic, bsp_blending_entire;
extern xge_widget     *sw_bsp_bl_nharmonic;
extern xge_int_widget bsp_blG1_minu, bsp_blG1_maxu, bsp_blG1_minv, bsp_blG1_maxv, 
                      bsp_blG2_minu, bsp_blG2_maxu, bsp_blG2_minv, bsp_blG2_maxv;
extern int            bsp_bl_opt_range[4];
extern int            bsp_bl_nkn1, bsp_bl_nkn2, bsp_bl_maxit;
extern double         sl_bsp_bl_param;
extern xge_int_widget bsp_bl_qknots1, bsp_bl_qknots2, bsp_bl_maxiter;
extern xge_widget     *bsp_bl_optimizeG1, *bsp_bl_optimizeG2;
extern xge_widget     *bsp_type_button;

extern boolean        sw_bsp_sproduct_equator, sw_bsp_sproduct_meridian,
                      sw_bsp_sproduct_eqrational,
                      sw_bsp_sproduct_eqclosed,
                      sw_bsp_sproduct_equniform,
                      sw_bsp_sproduct_merrational,
                      sw_bsp_sproduct_merclosed,
                      sw_bsp_sproduct_meruniform;
extern xge_string_ed  bsp_sproduct_eqname_ed, bsp_sproduct_mername_ed;
extern xge_int_widget intw_bsp_sproduct_eqdeg, intw_bsp_sproduct_merdeg;
extern int            bsp_eqdegree, bsp_merdegree;
extern GO_BSplineCurve *bsp_sproduct_eqmer;

    /* B-spline meshes */
extern xge_string_ed  bsm_name_ed;
extern xge_int_widget intw_bsmdeg;
extern xge_int_widget intw_bsm_vert0, intw_bsm_vert1,
                      intw_bsm_edge0, intw_bsm_edge1,
                      intw_bsm_fac0, intw_bsm_fac1;
extern int            bsm_vertex_num0, bsm_vertex_num1,
                      bsm_edge_num0, bsm_edge_num1,
                      bsm_facet_num0, bsm_facet_num1;
extern xge_int_widget intw_bsmdeg1, intw_bsmdeg2;
extern int            bsm_degree1, bsm_degree2;
extern xge_int_widget intw_bsm_simplify_x, intw_bsm_simplify_y, intw_bsm_simplify_z;
extern int            bsm_simplify_xyz[4];
extern xge_int_widget intw_bsm_decimate_iter;
extern int            bsm_decimate_iter;

extern xge_widget     *bsm_density;
extern boolean        bsm_sw_subdivision;
extern boolean        bsm_sw_blending;
extern boolean        bsm_sw_view_special;

extern boolean        bsm_sw_constr, bsm_sw_log_it, bsm_sw_shape_only,
                      bsm_sw_alt_multilevel;
extern double         sl_bsm_bl_param;
extern int            bsm_bl_nkn1, bsm_bl_nkn2, bsm_bl_maxit,
                      bsm_bl_nlev, bsm_bl_nbl;
extern xge_int_widget bsm_bl_qknots1, bsm_bl_qknots2, bsm_bl_maxiter,
                      bsm_bl_nlevels, bsm_bl_nblocks;
extern xge_widget     *bsm_bl_optimize, *bsc_mc_optimize;
extern boolean        bsm_sw_use_coarse;
extern xge_string_ed  bsm_coarse_editor;
extern xge_int_widget bsmnpthreads;
extern int            bsm_npthreads;

extern boolean        bsm_sw_data_add, bsm_sw_data_replace;

extern char           logfilename[];

    /* B-spline holes */

    /* colour specimen widget */
extern xge_widget     *colour_box;
extern double         colour_rgb[3];

    /* ray traced picture data */
extern boolean rendered_picture;
    /* lights & special stuff editing mode */
extern boolean picture_lights, editing_shapefunc, editing_lights,
               editing_pretrans, editing_camera;

    /* popup menu 14 - pre-transformation */
extern xge_string_ed pretrans_ed[12];
extern char          pretrans_str[12][MAX_PARAM_LGT+1];

xge_widget *InitTop00Menu ( xge_widget *prev );
int Top00MenuCallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *InitSide00Menu ( xge_widget *prev );
boolean ChangeSide00MenuWidth ( short h );
void SetupWin0StatusLine ( void );
void SetupWin0CommandLine ( void );
void SelectGeomTool ( char tool, boolean *sw );
void SetToolSwitches ( char tool, char coord );
int Side00MenuCallBack ( xge_widget *er, int msg, int key, short x, short y );

void InitXRenderer ( void );
void DestroyXRenderer ( void );
void SetRendererSwitches ( void );

void InitSide00bMenu ( void );
boolean InitRendering ( void );
void StartRendering ( void );
void StopRendering ( void );
void ContinueRendering ( void );
void SetupCameraParamWidgets ( void );
void UpdateCameraPosition ( CameraRecd *CPos, char which, double ang );
int Side00bMenuCallBack ( xge_widget *er, int msg, int key, short x, short y );

void InitSide00eMenu ( void );
int Side00eMenuCallBack ( xge_widget *er, int msg, int key, short x, short y );
void EnterRefLine ( point3d *p );

boolean VerifyNumParam ( char *txt, double *param );
boolean VerifyAngleParam ( char *txt, double *param );
void EscapeNumParam ( char *txt, double param );
void EscapeAngleParam ( char *txt, double param );

void BottomDisplayPoint2D ( int win, int pn, point2d *p, boolean onscreen );
void BottomDisplayPoint3D ( int win, int pn, point3d *p, boolean onscreen );
void BottomDisplayPoint4D ( int win, int pn, point4d *p, boolean onscreen );
void BottomDisplayPoint ( int win, int spdimen, int cpdimen,
                          int pn, double *pc, boolean onscreen );

xge_widget *InitBottom00Menu ( xge_widget *prev );
void ExecuteCommand0 ( void );
int Bottom00MenuCallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *InitGeom00Menu ( xge_widget *prev );
void DrawG00win3Dpar ( xge_widget *er, boolean onscreen );
void DrawG00win3Dpersp ( xge_widget *er, boolean onscreen );
void Geom00Win3DShowTransformation ( xge_3Dwind *ww );
void Geom00Win3DShowTransOrigin ( xge_3Dwind *ww );
void Geom00Win3DShowCameraPos ( CameraRecd *CPos );
void Geom00WinRestartRendering ( void );
int Geom00Win3DCallBack ( xge_widget *er, int msg, int key, short x, short y );
void DrawG00win2D ( xge_widget *er, boolean onscreen );
void Geom00Win2DShowTransformation ( xge_2Dwind *ww );
void Geom00Win2DShowTransOrigin ( xge_2Dwind *ww );
int Geom00Win2DCallBack ( xge_widget *er, int msg, int key, short x, short y );
void SetWin002D ( void );
void SetWin003D ( void );
int Geom00MenuCallBack ( xge_widget *er, int msg, int key, short x, short y );

void SetCurrentDir ( void );
xge_widget *InitPopup00 ( void );
void PreparePopup01 ( void );
void PreparePopup02 ( void );
boolean FilenameCorrect ( char *fn );
int Popup00CallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *InitPopup01 ( void );
void Popup01ChangeDir ( void );
void Popup01ChangeDirAlt ( short x );
void Popup01CameraReader ( void *usrdata, int ident, CameraRecd *camera );
void Popup01OpenFile ( void );
int Popup01CallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *InitPopup02 ( void );
void Popup02ChangeDir ( void );
void Popup02ChangeDirAlt ( short x );
boolean WriteOtherData ( void *usrdata );
void Popup02SaveFile ( void );
int Popup02CallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *InitPopup03 ( void );
int Popup03CallBack ( xge_widget *er, int msg, int key, short x, short y );


xge_widget *InitTop10Menu ( xge_widget *prev );
int Top10MenuCallBack ( xge_widget *er, int msg, int key, short x, short y );

void InitNameEditor ( xge_string_ed *ed, char *name );
xge_widget *InitSide10Menu ( xge_widget *prev );
boolean ChangeSide10MenuWidth ( short h );
void SetupWin1StatusLine ( void );
void SetupWin1CommandLine ( void );
int Side10MenuCallBack ( xge_widget *er, int msg, int key, short x, short y );

void InitSide10Menu_Bezc ( void );
boolean ChangeSide10MenuWidth_Bezc ( short h );
void SetupBezierCurveWidgets ( GO_BezierCurve *obj );
int Side10MenuBezcCallBack ( xge_widget *er, int msg, int key, short x, short y );

void InitSide10Menu_Bezp ( void );
boolean ChangeSide10MenuWidth_Bezp ( short h );
void SetupBezierPatchWidgets ( GO_BezierPatch *obj );
int Side10MenuBezpCallBack ( xge_widget *er, int msg, int key, short x, short y );

void InitSide10Menu_BSc ( void );
boolean ChangeSide10MenuWidth_BSc ( short h );
void SetupBSplineCurveWidgets ( GO_BSplineCurve *obj );
int Side10MenuBscCallBack ( xge_widget *er, int msg, int key, short x, short y );
boolean MengerCurvatureOptimizationPrepareData ( GO_BSplineCurve *obj );
void InitMengerCurvatureOptimization ( void );

void InitSide10Menu_BSp ( void );
boolean ChangeSide10MenuWidth_BSp ( short h );
void SetupBSplinePatchWidgets ( GO_BSplinePatch *obj );
void BSPatchSetupGeneralOptions ( GO_BSplinePatch *obj );
void BSPatchSetupSphericalOptions ( GO_BSplinePatch *obj );
void BSPatchSetupBlG1Options ( GO_BSplinePatch *obj );
void BSPatchSetupBlG2Options ( GO_BSplinePatch *obj );
GO_BSplineMesh *BlendingMeshOptimizationFindCoarse ( GO_BSplineMesh *obj );
boolean BlendingPatchOptimizationPrepareData ( GO_BSplinePatch *obj );
void InitBlendingPatchOptimization ( void );
void SuggestOptLevels ( GO_BSplineMesh *obj, boolean multilevel );
int Side10MenuBspCallBack ( xge_widget *er, int msg, int key, short x, short y );

void InitSide10Menu_BSm ( void );
boolean ChangeSide10MenuWidth_BSm ( short h );
void SetupBSplineMeshVEFnum ( GO_BSplineMesh *obj );
void SetupBSplineMeshWidgets ( GO_BSplineMesh *obj );
boolean BlendingMeshOptimizationPrepareData ( GO_BSplineMesh *obj );
void InitBlendingMeshOptimization ( void );
void BlendingMeshOptSpecialPatches ( void );
void Side10MenuBsmResize ( short x, short y );
int Side10MenuBsmCallBack ( xge_widget *er, int msg, int key, short x, short y );

void InitSide10Menu_BSh ( void );
boolean ChangeSide10MenuWidth_BSh ( short h );
void SetupBSplineHoleWidgets ( GO_BSplineHole *obj );
int Side10MenuBshCallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *InitBottom10Menu ( xge_widget *prev ); 
void ExecuteCommand1 ( void );
int Bottom10MenuCallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *InitGeom10Menu ( xge_widget *prev );
void Geom10MenuResize ( void );
int Geom10MenuCallBack ( xge_widget *er, int msg, int key, short x, short y );
void DrawG10win2D ( xge_widget *er, boolean onscreen );
int Geom10win2DCallBack ( xge_widget *er, int msg, int key, short x, short y );
void Geom10winKNSetKnots ( xge_KnotWind *kwin, int degree, int lastknot,
                           double *knots, boolean closed );
int Geom10winKNCallBack ( xge_widget *er, int msg, int key, short x, short y );
void DrawG10winT2KNKnots ( xge_T2KnotWind *kwin );
void DrawG10winT2KNKnotLines ( xge_T2KnotWind *kwin );
void DrawG10winT2KNDomain ( xge_T2KnotWind *kwin );
void DrawG10winT2KNDomainNet ( xge_T2KnotWind *kwin );
void DrawG10winTrimmedDomain ( xge_T2KnotWind *kwin );
void DrawG10winT2KN ( xge_widget *er, boolean onscreen );
void Geom10winT2SetKnots ( xge_T2KnotWind *kwin,
                           int degreeu, int lastknotu, double *knotsu, boolean closedu,
                           int degreev, int lastknotv, double *knotsv, boolean closedv );
int Geom10winT2KNCallBack ( xge_widget *er, int msg, int key, short x, short y );
void SetGeomWin10Empty ( void );
void SetGeomWin10Knotw ( GO_BSplineCurve *go );
void SetGeomWin10T2Knotw ( GO_BSplinePatch *go );
void SetGeomWin10Win2D ( GO_BSplineHole *go );
void SetGeomWin10SPrEqMer ( GO_BSplineCurve *eqmer );
void DrawG10win2Deqmer ( xge_widget *er, boolean onscreen );
int Geom10win2DeqmerCallBack ( xge_widget *er, int msg, int key, short x, short y );
int Geom10winKNeqmerCallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *InitPopup10 ( void );
int Popup10CallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *InitPopup11 ( void );
void SetupPopup11Bezc ( void );
void SetupPopup11Bezp ( void );
void SetupPopup11Bsc ( void );
void SetupPopup11Bsp ( void );
void SetupPopup11Bsm ( void );
void SetupPopup11Bsh ( void );
void AddANewObject ( int obj_type );
int Popup11CallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *InitPopup12 ( void );
boolean SetupObjectNameList ( void );
void DeleteObjectNameList ( void );
void UpdateObjectNameList ( void );
void CleanupPopup12 ( void );
int Popup12CallBack ( xge_widget *er, int msg, int key, short x, short y );

void ColourSpecimenRedraw ( xge_widget *er, boolean onscreen );
xge_widget *InitPopup13 ( void );
int Popup13CallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *InitPopup14 ( void );
void GetPreTransformation ( void );
void CleanupPopup14 ( void );
int Popup14CallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *InitPopup15 ( void );
int Popup15CallBack ( xge_widget *er, int msg, int key, short x, short y );

void SetupObjectSpecificMenusWin1 ( geom_object *obj );
void SetupObjectSpecificMenus ( geom_object *obj );
void SetObjectEditMenu ( geom_object *obj );
void SetObjectViewMenu ( geom_object *obj );
void SetObjectDataMenu ( geom_object *obj );
void SetObjectOptionsMenu ( geom_object *obj );

void SetStatusText ( char *text, boolean onscreen );
void NotifyParam1 ( double param );
void NotifyParam2 ( double param );

void RedrawGeom00Win ( void );
void SetCircleCursor ( int win );
void OpenPopup ( xge_widget *er, boolean allwin );
void InitWindow0Widgets ( void );
void InitWindow1Widgets ( void );
boolean ProcessCMDLineParameters ( int argc, char *argv[] );
void init_program ( int argc, char *argv[] );
void destroy_program ( void );

void ResizeWindow0 ( short x, short y );
int Window0CallBack ( int msg, int key, short x, short y );

void ResizeWindow1 ( short x, short y );
int Window1CallBack ( int msg, int key, short x, short y );

void HandleChildFailure ( void );
int ProcessIdleCommand ( int key, short x, short y );
int NonWidgetCallBack ( int msg, int key, short x, short y );
int CallBack ( xge_widget *er, int msg, int key, short x, short y );

