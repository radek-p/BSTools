
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define IDLE_COMMAND_START_RENDERING  1
#define IDLE_COMMAND_RENDER           2
#define IDLE_COMMAND_SEND_OPTIONS     3
#define IDLE_COMMAND_SEND_PATCH       4
#define IDLE_COMMAND_SEND_CONSTRAINTS 5
#define IDLE_COMMAND_RECEIVE_PATCH    6

#define MAX_PATH_LGT     1024
#define MAX_FILENAME_LGT   64

#define STATUS_LINE_LENGTH ((xge_MAX_WIDTH)/6)

extern int  program_argc;
extern char **program_argv;

extern int win0, win1;
extern xge_widget *menu0, *menu1;
extern xge_widget *menu00list, *menu01list, *menu02list, *menu03list;
extern xge_widget *popup00, *popup01, *popup02, *popup03, *popup04;
extern xge_int_widget ddeg_u, ddeg_v, ddens_u, ddens_v;

extern xge_widget *domwind, *menu2, *menu3;
extern xge_widget *menu4, *menu5, *menu6, *menu7, *menu8, *menu9,
                  *menu10, *menu11, *menu12, *menu13, *menu14;

extern xge_widget *menu60list, *menu61list, *menu62list;
extern xge_widget *popup10, *popup11, *popup12;

extern xge_listbox filelist, dirlist;
extern xge_string_ed filename_editor;

extern char       *InfoMsg[5];

extern const char file_filter[];
extern const char file_ext[];
extern char current_directory[];
extern char filename[];

extern xge_widget *status0, *status0sw;
extern boolean status_0_sw;
extern char statustext0[];

extern xge_widget *renderbtn0, *renderbtn1;

extern boolean edit_light_sw[4];
extern boolean edit_reflection_frame;
extern boolean edit_highlight_frame;
extern boolean edit_sections_frame;

extern boolean swShadows, swAntialias;
extern boolean swind_picture;

extern double dfs_factor;
extern double cf_range[2];

/* light directions and intensities - the last is the ambient intensity */
extern vector3d light_dir[R_NLIGHTS];
extern double   light_int[R_NLIGHTS+1];

extern xge_widget *status1, *status1sw;
extern boolean status_1_sw;
extern char statustext1[];

/* iteration number limit for nonlinear blending surfaces */
extern xge_int_widget wdg_blending_lmt_iter;
extern xge_widget *bl_optbtn;
extern int blending_quad1, blending_quad2;
extern xge_int_widget wdg_blending_quad1, wdg_blending_quad2;

extern int ipc_state;
extern boolean constraints_sent;

extern xge_widget *bl_trihsw;
extern xge_int_widget blending_opt_range[4];
extern xge_string_ed blending_trans_editor[9];
extern char blending_trans_str[9][13];
extern boolean sw_show_steps;

/* ///////////////////////////////////////////////////////////////////////// */
extern boolean child_ready;

/* ///////////////////////////////////////////////////////////////////////// */
void SetBox2d ( Box2d *box, double x0, double x1, double y0, double y1 );
void SetBox3d ( Box3d *box,       
                double x0, double x1, double y0, double y1, double z0, double z1 );

void SetupDomainMapping ( void );
void ZoomDomainWindow ( double scf );
boolean PanDomainWindow ( short x, short y );
void ResetDomainMapping ( void );
void FindDomainMapping ( void );
void ResizeDomainWinContents ( short w, short h, short x, short y );

void RedrawGeomWindows ( boolean w0, boolean w1 );
void RysujOkno ( xge_widget *er, boolean onscreen );
void RysujPOkno ( xge_widget *er, boolean onscreen );
void RysujDOkno ( xge_widget *er, boolean onscreen );
void RysujSPROkno ( xge_widget *er, boolean onscreen );

void init_edwin ( void );
void destroy_edwin ( void );

void BreakSprBind ( void );
void BreakNLBlending ( void );
void ProcessIdleCommand ( int key, short x, short y );
void ProcessChildMessage ( int msg, int size );
int ProcessOtherMsg ( int msg, int key, short x, short y );
void OpenFile ( void );
void UpdateBlendingRangeWidgets ( void );
int Win0CallBack ( xge_widget *er, int msg, int key, short x, short y );
int Win1CallBack ( xge_widget *er, int msg, int key, short x, short y );
int CallBack ( xge_widget *er, int msg, int key, short x, short y );

boolean FilenameCorrect ( char *filename );
boolean WritePatchAttributes ( void *usrdata );
boolean SaveBSPatch ( char *filename );
void BSPatchReader ( void *userData, const char *name, int ident,
                     int udeg, int lastknotu, const double *knotsu, 
                     int vdeg, int lastknotv, const double *knotsv, 
                     boolean closed_u, boolean closed_v, 
                     int pitch, const point4d *cpoints, 
                     int spdimen, boolean rational );
void CPMarkReader ( void *userData, int ncp, unsigned int *mk );
boolean ReadBSPatch ( char *filename );

void ResizeWinStatus ( int win );
void StatusLineOnOff ( int win );
void SetStatusLineText ( int win, const char *text, boolean onscreen );
void Notify2DTrans ( int win, xge_2Dwind *_2Dwin );
void Notify2DTransChange ( int win, xge_2Dwind *_2Dwin );
void Notify3DTrans ( int win, xge_3Dwind *_3Dwin );
void Notify3DTransChange ( int win, xge_3Dwind *_3Dwin );
void NotifyDoubleNumber ( int win, double x, boolean onscreen );

void StartRendering ( void );
void ContRendering ( void );
void BreakRendering ( boolean redraw );

void SetKWindNKnots ( void );

void ResetBlendingOptTrans ( void );
boolean EnterBlendingOptTransCoefficient ( int id );

boolean MakeTheChildProcess ( void );

