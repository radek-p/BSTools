
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#define IDLE_COMMAND_FIND_SURFACE1   1
#define IDLE_COMMAND_FIND_SURFACE2   2
#define IDLE_COMMAND_START_RENDERING 3
#define IDLE_COMMAND_RENDER          4

#define MAX_PATH_LGT       1024
#define MAX_FILENAME_LGT   64

#define STATUS_LINE_LENGTH ((xge_MAX_WIDTH-110)/6)

/* the #definitions below must match these in edwidgets.h */
#define MAX_CONSTR_SWITCH_ROWS     16
#define MAX_CONSTR_SWITCH_COLS     16
#define NUM_CONSTR_SWITCHES       257  /* 16*16+1 */


extern int win0, win1;
extern xge_widget *menu0, *menu1, *menu2;
extern xge_widget *menu01alist, *menu01blist, *menu01clist,
                  *menu01dlist, *menu01elist, *menu01flist;
extern xge_widget *popup00, *popup01, *popup02, *popup03, *popup04;
extern xge_widget *status0, *status0sw, *renderbtn0, *renderbtn1;
extern char statustext0[];
extern boolean status_0_sw;
extern xge_int_widget hole_k_intw;

/* constraint switches - to be reconfigured whenever necessary */
extern xge_widget *menu01cscrolled;
extern xge_scroll_widget scroll_constr_sw;
extern xge_widget *constr_sw[NUM_CONSTR_SWITCHES];
extern boolean constraint[NUM_CONSTR_SWITCHES];
extern xge_widget *constr_row_text[MAX_CONSTR_SWITCH_ROWS];
extern boolean constraints1, constraints2,
               constraints1st, constraints2nd, constraintsZero;

extern xge_listbox   filelist, dirlist;
extern xge_string_ed filename_editor;
extern char current_directory[];
extern char filename[];
extern const char file_ext[];
extern const char file_filter[];

extern xge_widget *menu3, *menu4, *menu5;
extern xge_widget *menu14alist, *menu14a1list, *menu14a2list,
                  *menu14blist, *menu14clist, *menu14dlist;
extern xge_widget *knotwind;
extern xge_widget *status1, *status1sw;
extern char statustext1[];
extern boolean status_1_sw;

extern boolean sw_opt_1, sw_opt_2;
extern xge_int_widget opt_order_1, opt_nk_1, opt_m1_1, opt_m2_1;
extern xge_int_widget opt_order_2, opt_nk_2, opt_m1_2, opt_m2_2;

extern boolean view_constraints_frame;
extern boolean view_surf_cp;
extern boolean view_surf_spatches;
extern boolean view_surf_1, view_surf_2;
extern boolean view_surf_numbers;

extern boolean view_dom_cp;
extern boolean view_dom_spatches;
extern boolean view_dom_patches1, view_dom_patches2;
extern boolean view_dom_numbers;

extern boolean edit_reflection_frame;
extern boolean edit_highlight_frame;
extern boolean edit_sections_frame;
extern boolean swShadows, swAntialias;

/* light directions and intensities - the last is the ambient intensity */
extern boolean edit_light_sw[R_NLIGHTS];
extern vector3d light_dir[R_NLIGHTS];
extern double   light_int[R_NLIGHTS+1];

void RysujOkno ( xge_widget *er, boolean onscreen );
void RysujPOkno ( xge_widget *er, boolean onscreen );

void RysujDomOkno ( xge_widget *er, boolean onscreen );


typedef struct {
          xge_widget *er;
          double     akm, bkm, umin, umax, knotscf, knotshf;
          double     *knots;
          char       hole_k;
          short      xx;
          boolean    panning, display_coord;
          char       current_seq, current_knot;
        } ghKnotWind;

void DrawKnotWind ( xge_widget *er, boolean onscreen );
void DrawKnotWindCursorPos ( ghKnotWind *knw );
boolean KnotWindMsg ( xge_widget *er, int msg, int key, short x, short y );
void KnotWindInitMapping ( ghKnotWind *knw, double umin, double umax );
void KnotWindSetHoleK ( ghKnotWind *knw, char hole_k, double *knots );
xge_widget *NewGHKnotWind ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            char hole_k, ghKnotWind *knw, double *knots );
short MapKnot ( ghKnotWind *knw, double x );
double UnmapKnot ( ghKnotWind *knw, short xi );
boolean FindNearestKnot ( ghKnotWind *knw, short x, short y );
boolean SetKnot ( ghKnotWind *knw, short x );
void KnotWindPan ( ghKnotWind *knw, short dxi );
void KnotWindZoom ( ghKnotWind *knw, double scf );


void init_edwin ( void );
void destroy_edwin ( void );
void SetupWindow0Widgets ( void );
void SetupWindow1Widgets ( void );
boolean ConstructConstraintSwitches ( void );
void ConfigureConstraintWidgets ( boolean reset );
void SwitchTheConstraint ( int constrno );

int ProcessKey ( int key );
void ProcessIdleCommand ( int key, short x, short y );
boolean SetupHoleSides ( int key );
void UpdateSurfPicture ( void );
int Win0CallBack ( xge_widget *er, int msg, int key, short x, short y );
int Win1CallBack ( xge_widget *er, int msg, int key, short x, short y );
int CallBack ( xge_widget *er, int msg, int key, short x, short y );

boolean FilenameCorrect ( const char *fn );
void OpenFile ( const char *fn );
void SaveFile ( const char *fn );

void ResizeWinStatus ( int win );
void SetStatusLineText ( int win, const char *text, boolean onscreen );
void StatusLineOnOff ( int win );
void Notify2DTrans ( int win, xge_2Dwind *_2Dwin );
void Notify2DTransChange ( int win, xge_2Dwind *_2Dwin );
void Notify3DTrans ( int win, xge_3Dwind *_3Dwin );
void Notify3DTransChange ( int win, xge_3Dwind *_3Dwin );
void NotifyDoubleNumber ( int win, double x, boolean onscreen );

void StartRendering ( void );
void ContRendering ( void );
void BreakRendering ( boolean redraw );

void DisplayBasisInfo ( void );

