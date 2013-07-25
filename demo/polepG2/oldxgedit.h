
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_WINDOWS    8

#define MAX_WIDTH   1024
#define MAX_HEIGHT   768
#define WIDTH        480
#define HEIGHT       360

#define msg_NULL             0
#define msg_ENTERING     0x100
#define msg_EXITING      0x101
#define msg_INIT         0x102
#define msg_KEY          0x103
#define msg_MMOVE        0x104
#define msg_MCLICK       0x105

#define mouse_LBUTTON_DOWN    (1 << 0)
#define mouse_LBUTTON_CHANGE  (1 << 1)
#define mouse_RBUTTON_DOWN    (1 << 2)
#define mouse_RBUTTON_CHANGE  (1 << 3)
#define mouse_MBUTTON_DOWN    (1 << 4)
#define mouse_MBUTTON_CHANGE  (1 << 5)

#define BORDER_WIDTH 2

#define RECT_NONE -1

extern int winnum;
extern unsigned int mouse_buttons;
extern int  mouse_x, mouse_y, xx, yy;

extern Display    *thedisplay;
extern Window     thewindow;
extern Pixmap     thepixmap;
extern int        thescreen;
extern GC         thegc;
extern Visual     *thevisual;
extern XSizeHints thehints;
extern Cursor     thecursor0, thecursor1;
extern KeySym     thekeysym;
extern XEvent     theevent;

extern unsigned int current_width, current_height;  /* window dimensions */

extern char       *p_name;
extern char       done;
extern int        prevx, prevy;
extern int        slidebarid;

extern Colormap   thecolormap;
extern int        ncolors, nplanes;
extern unsigned long foreground, background;
extern unsigned long palette[];   /* my private palette */

#define c_black      (palette[0x00])
#define c_dk_red     (palette[0x04])
#define c_red        (palette[0x06])
#define c_lt_red     (palette[0x07])
#define c_dk_green   (palette[0x10])
#define c_green      (palette[0x28])
#define c_lt_green   (palette[0x30])
#define c_llt_green  (palette[0x79])
#define c_dk_yellow  (palette[0x36])
#define c_yellow     (palette[0x3F])
#define c_dk_blue    (palette[0x40])
#define c_dk_magenta (palette[0x43])
#define c_lt_magenta (palette[0xD7])
#define c_dk_grey    (palette[0x52])
#define c_lt_yellow  (palette[0x7F])
#define c_blue       (palette[0x80])
#define c_lt_blue    (palette[0xA2])
#define c_cyan       (palette[0xF8])
#define c_lt_grey    (palette[0xAD])
#define c_white      (palette[0xFF])


typedef struct ed_rect {
          int id;
          int w, h, x, y;
          void (*msgproc) ( struct ed_rect*, int, int, int, int );
          void (*redraw) ( struct ed_rect* );
        } ed_rect;
typedef void (*message_proc) ( ed_rect *er, int msg, int key, int x, int y );
typedef void (*redraw_proc) ( ed_rect *er );


/* States of the program. Other may be defined in applications. */

#define STATE_NOTHING     0
#define STATE_MOVINGSLIDE 1
#define STATE_MESSAGE     2

extern int     stan;
extern ed_rect *focus;
extern boolean notinfocus;


/* xgedit.c procedures */

int     NewWindow ( char *p_name );
boolean SetWindow ( int win );
int     CurrentWindow ( void );

void redraw ( void );
void redraw_all ( void );
void DrawEmpty ( ed_rect *er );
void EmptyMsg ( ed_rect *er, int msg, int key, int x, int y );
void DrawMenu ( ed_rect *er );
void MenuMsg ( ed_rect *er, int msg, int key, int x, int y );
void DrawSwitch ( ed_rect *er );
void SwitchMsg ( ed_rect *er, int msg, int key, int x, int y );
void DrawButton ( ed_rect *er );
void ButtonMsg ( ed_rect *er, int msg, int key, int x, int y );
void DrawSlidebar ( ed_rect *er );
void SlidebarMsg ( ed_rect *er, int msg, int key, int x, int y );
void DrawText ( ed_rect *er );

float LinSlidebarValue ( float xmin, float xmax, float t );
float LogSlidebarValue ( float xmin, float xmax, float t );

void DisplayErrorMessage ( char *message );
void DisplayInfoMessage ( char **msglines );

void BoundPoint ( ed_rect *er, int *x, int *y );
boolean RectanglesIntersect ( int wa, int ha, int xa, int ya,
                              int wb, int hb, int xb, int yb );

void SetupEdRect ( ed_rect *edr, int en, int id, int w, int h, int x, int y,
                   message_proc msgproc, redraw_proc redraw );
void SetWinEdRect ( int rect_num, ed_rect *edr );

void MakeTCPalette ( void );
void GetWindowSize ( void );
void init ( int argc, char *argv[] );
void cleanup ( void );

boolean PointInRect ( ed_rect *edr, int x, int y );

void dispatch_message ( unsigned int msg, unsigned int key, int x, int y );
void get_message ( unsigned int *msg, unsigned int *key, int *x, int *y );


/* ////////////////////////////////////////////////////////////////////////// */
/* user supplied routines, outside the file xgedit.c: */

void init_edwin ( int argc, char *argv[] );
void destroy_edwin ( void );
void resize_edwin ( void );
void process_key ( unsigned int *key );
void GetMenuData ( ed_rect *er, int *menu_rect_num, ed_rect **menu );
void GetSwitchData ( ed_rect *er,
                     char **title, boolean **switchvar,
                     void (**switchproc)(ed_rect *er) );
void GetButtonData ( ed_rect *er,
                     char **title, void (**buttonproc)(ed_rect *er) );
void GetSlidebarData ( ed_rect *er,
                       float **slipos, void (**slideproc)(ed_rect *er) );
void GetTextTitle ( ed_rect *er,
                    char **title );
void ProcessOtherEvent ( void );

