
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef XGEDIT_H
#define XGEDIT_H

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#define USE_XEXT_SHAPE
#ifdef USE_XEXT_SHAPE
#include <X11/extensions/shape.h>
#endif

#ifndef _LIBC_LIMITS_H_
#include <limits.h>
#endif

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif
#ifndef PKNUM_H
#include "pknum.h"
#endif
#ifndef PKGEOM_H
#include "pkgeom.h"
#endif
#ifndef CAMERA_H
#include "camera.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* sometimes automatic screen aspect detection produces wrong results; */
/* that being the case, comment out the #define XGE_AUTO_ASPECT line,  */
/* perhaps a better solution will some day be found (dealing also      */
/* with virtual displays of size other than physical screens)          */
#define XGE_AUTO_ASPECT
#ifndef XGE_AUTO_ASPECT
#define XGE_DEFAULT_ASPECT     1.0
#endif

#define xge_MAX_WINDOWS          8
#define xge_MAX_CURSORS         16

/* widget input focus stack capacity */
#define xge_FOCUS_DEPTH          8

/* maximal and initial window sizes */
#define xge_MAX_WIDTH         1280
#define xge_MAX_HEIGHT         960
#define xge_WIDTH              480
#define xge_HEIGHT             360
/* maximal special widget size */
#define xge_MAX_SPECIAL_WIDTH  360
#define xge_MAX_SPECIAL_HEIGHT 360

/* size of characters in the fixed font */
#define xge_CHAR_WIDTH           6
#define xge_CHAR_HEIGHT         13

/* tolerance of pointing for some widgets */
#define xge_RECT_NONE           -1
#define xge_MINDIST              8

/* limit for text editing widgets */
#define xge_MAX_STRING_LENGTH  512

/* ///////////////////////////////////////////////////////////////////////// */
/* mouse button states */
#define xgemouse_LBUTTON_DOWN    (1 << 0)
#define xgemouse_LBUTTON_CHANGE  (1 << 1)
#define xgemouse_RBUTTON_DOWN    (1 << 2)
#define xgemouse_RBUTTON_CHANGE  (1 << 3)
#define xgemouse_MBUTTON_DOWN    (1 << 4)
#define xgemouse_MBUTTON_CHANGE  (1 << 5)
#define xgemouse_WHEELFW_DOWN    (1 << 6)
#define xgemouse_WHEELFW_CHANGE  (1 << 7)
#define xgemouse_WHEELBK_DOWN    (1 << 8)
#define xgemouse_WHEELBK_CHANGE  (1 << 9)

/* ///////////////////////////////////////////////////////////////////////// */
/* message codes */
#define xgemsg_NULL                              0
#define xgemsg_INIT                          0x100

/* user action messages */
#define xgemsg_KEY                           0x101
#define xgemsg_SPECIAL_KEY                   0x102
#define xgemsg_MMOVE                         0x103
#define xgemsg_MCLICK                        0x104
#define xgemsg_OTHEREVENT                    0x105

/* commands sent by libxgedit procedures */
#define xgemsg_ENTERING                      0x106
#define xgemsg_EXITING                       0x107
#define xgemsg_RESIZE                        0x108
#define xgemsg_MOVE                          0x109

#define xgemsg_BUTTON_COMMAND                0x10A
#define xgemsg_SWITCH_COMMAND                0x10B
#define xgemsg_SLIDEBAR_COMMAND              0x10C
#define xgemsg_SLIDEBAR2_COMMAND             0x10D
#define xgemsg_DIAL_COMMAND                  0x10E
#define xgemsg_TEXT_EDIT_VERIFY              0x10F
#define xgemsg_TEXT_EDIT_ENTER               0x110
#define xgemsg_TEXT_EDIT_ESCAPE              0x111
#define xgemsg_INT_WIDGET_COMMAND            0x112
#define xgemsg_LISTBOX_ITEM_SET              0x113
#define xgemsg_LISTBOX_ITEM_PICK             0x114
#define xgemsg_LISTBOX_EXCHANGE              0x115
#define xgemsg_QUATROTBALL_COMMAND           0x116
#define xgemsg_TEXT_WIDGET_CLICK             0x117

#define xgemsg_2DWIN_RESIZE                  0x118
#define xgemsg_2DWIN_PROJCHANGE              0x119
#define xgemsg_2DWIN_PICK_POINT              0x11A
#define xgemsg_2DWIN_MOVE_POINT              0x11B
#define xgemsg_2DWIN_SELECT_POINTS           0x11C
#define xgemsg_2DWIN_UNSELECT_POINTS         0x11D
#define xgemsg_2DWIN_SPECIAL_SELECT          0x11E
#define xgemsg_2DWIN_SPECIAL_UNSELECT        0x11F
#define xgemsg_2DWIN_CHANGE_TRANS            0x120
#define xgemsg_2DWIN_SAVE_POINTS             0x121
#define xgemsg_2DWIN_TRANSFORM_POINTS        0x122
#define xgemsg_2DWIN_TRANSFORM_SPECIAL       0x123
#define xgemsg_2DWIN_FIND_REFBBOX            0x124
#define xgemsg_2DWIN_UNDO                    0x125
#define xgemsg_2DWIN_KEY                     0x126
#define xgemsg_2DWIN_ERROR                   0x127

#define xgemsg_3DWIN_RESIZE                  0x128
#define xgemsg_3DWIN_PROJCHANGE              0x129
#define xgemsg_3DWIN_PICK_POINT              0x12A
#define xgemsg_3DWIN_MOVE_POINT              0x12B
#define xgemsg_3DWIN_SELECT_POINTS           0x12C
#define xgemsg_3DWIN_UNSELECT_POINTS         0x12D
#define xgemsg_3DWIN_SPECIAL_SELECT          0x12E
#define xgemsg_3DWIN_SPECIAL_UNSELECT        0x12F
#define xgemsg_3DWIN_CHANGE_TRANS            0x130
#define xgemsg_3DWIN_SAVE_POINTS             0x131
#define xgemsg_3DWIN_TRANSFORM_POINTS        0x132
#define xgemsg_3DWIN_TRANSFORM_SPECIAL       0x133
#define xgemsg_3DWIN_FIND_REFBBOX            0x134
#define xgemsg_3DWIN_UNDO                    0x135
#define xgemsg_3DWIN_KEY                     0x136
#define xgemsg_3DWIN_ERROR                   0x137

#define xgemsg_KNOTWIN_CHANGE_KNOT           0x138
#define xgemsg_KNOTWIN_INSERT_KNOT           0x139
#define xgemsg_KNOTWIN_REMOVE_KNOT           0x13A
#define xgemsg_KNOTWIN_CHANGE_ALTKNOT        0x13B
#define xgemsg_KNOTWIN_INSERT_ALTKNOT        0x13C
#define xgemsg_KNOTWIN_REMOVE_ALTKNOT        0x13D
#define xgemsg_KNOTWIN_MCLICK                0x13E
#define xgemsg_KNOTWIN_MMOVE                 0x13F
#define xgemsg_KNOTWIN_CHANGE_MAPPING        0x140
#define xgemsg_KNOTWIN_ERROR                 0x141

#define xgemsg_T2KNOTWIN_RESIZE              0x142
#define xgemsg_T2KNOTWIN_PROJCHANGE          0x143
#define xgemsg_T2KNOTWIN_CHANGE_KNOT_U       0x144
#define xgemsg_T2KNOTWIN_CHANGE_KNOT_V       0x145
#define xgemsg_T2KNOTWIN_INSERT_KNOT_U       0x146
#define xgemsg_T2KNOTWIN_INSERT_KNOT_V       0x147
#define xgemsg_T2KNOTWIN_REMOVE_KNOT_U       0x148
#define xgemsg_T2KNOTWIN_REMOVE_KNOT_V       0x149
#define xgemsg_T2KNOTWIN_CHANGE_ALTKNOT_U    0x14A
#define xgemsg_T2KNOTWIN_CHANGE_ALTKNOT_V    0x14B
#define xgemsg_T2KNOTWIN_INSERT_ALTKNOT_U    0x14C
#define xgemsg_T2KNOTWIN_INSERT_ALTKNOT_V    0x14D
#define xgemsg_T2KNOTWIN_REMOVE_ALTKNOT_U    0x14E
#define xgemsg_T2KNOTWIN_REMOVE_ALTKNOT_V    0x14F
#define xgemsg_T2KNOTWIN_SELECT_POINTS       0x150
#define xgemsg_T2KNOTWIN_UNSELECT_POINTS     0x151
#define xgemsg_T2KNOTWIN_CHANGE_MAPPING      0x152
#define xgemsg_T2KNOTWIN_ERROR               0x153

/* message sent when popups are removed */
#define xgemsg_POPUP_REMOVED                 0x154
#define xgemsg_POPUPS_REMOVED                0x155

/* message sent after an error, warning or info message has been dismissed */
#define xgemsg_USER_MESSAGE_DISMISSED        0x156

/* message sent to the application by itself, to be delivered */
/* when there is nothing else to do */
#define xgemsg_IDLE_COMMAND                  0x157

/* message sent to the application by a child process */
#define xgemsg_CHILD_MESSAGE                 0x158
/* message sent after the child process has been terminated */
#define xgemsg_CHILD_FAILURE                 0x159

/* additional application messages must be greater than xgemsg_LAST_MESSAGE */
#define xgemsg_LAST_MESSAGE xgemsg_CHILD_FAILURE

/* ///////////////////////////////////////////////////////////////////////// */
/* States of the program. Others may be defined in applications. */
#define xgestate_NOTHING                      0
#define xgestate_BUTTON_DEFAULT               1
#define xgestate_BUTTON_COMBO_0               2
#define xgestate_BUTTON_COMBO_1               3
#define xgestate_BUTTON_INACTIVE              4
#define xgestate_MOVINGSLIDE                  5
#define xgestate_MOVINGSLIDE2A                6
#define xgestate_MOVINGSLIDE2B                7
#define xgestate_TURNINGDIAL                  8
#define xgestate_QUATROT_TURNING1             9
#define xgestate_QUATROT_TURNING2            10
#define xgestate_QUATROT_TURNING3            11
#define xgestate_MESSAGE                     12
#define xgestate_RESIZING_X                  13
#define xgestate_RESIZING_Y                  14
#define xgestate_RESIZING_XY                 15
#define xgestate_TEXT_EDITING                16
#define xgestate_2DWIN_MOVINGPOINT           17
#define xgestate_2DWIN_PANNING               18
#define xgestate_2DWIN_ZOOMING               19
#define xgestate_2DWIN_SELECTING             20
#define xgestate_2DWIN_UNSELECTING           21
#define xgestate_2DWIN_MOVING_GEOM_WIDGET    22
#define xgestate_2DWIN_USING_GEOM_WIDGET     23
#define xgestate_2DWIN_ALTUSING_GEOM_WIDGET  24
#define xgestate_2DWIN_USING_SPECIAL_WIDGET  25
#define xgestate_3DWIN_MOVINGPOINT           26
#define xgestate_3DWIN_PARPANNING            27
#define xgestate_3DWIN_PARZOOMING            28
#define xgestate_3DWIN_TURNING_VIEWER        29
#define xgestate_3DWIN_PANNING               30
#define xgestate_3DWIN_ZOOMING               31
#define xgestate_3DWIN_SELECTING             32
#define xgestate_3DWIN_UNSELECTING           33
#define xgestate_3DWIN_MOVING_GEOM_WIDGET    34
#define xgestate_3DWIN_USING_GEOM_WIDGET     35
#define xgestate_3DWIN_ALTUSING_GEOM_WIDGET  36
#define xgestate_3DWIN_USING_SPECIAL_WIDGET  37
#define xgestate_KNOTWIN_MOVINGKNOT          38
#define xgestate_KNOTWIN_PANNING             39
#define xgestate_KNOTWIN_ZOOMING             40
#define xgestate_T2KNOTWIN_MOVINGKNOT_U      41
#define xgestate_T2KNOTWIN_MOVINGKNOT_V      42
#define xgestate_T2KNOTWIN_MOVING_POINT      43
#define xgestate_T2KNOTWIN_PANNING           44
#define xgestate_T2KNOTWIN_ZOOMING           45
#define xgestate_T2KNOTWIN_SELECTING         46
#define xgestate_T2KNOTWIN_UNSELECTING       47
#define xgestate_LISTBOX_PICKING             48
/* additional application states must be greater than xge_LAST_STATE */
#define xgestate_LAST  xgestate_LISTBOX_PICKING

/* ///////////////////////////////////////////////////////////////////////// */
/* cursor identifiers */
#define xgeCURSOR_CROSSHAIR xgecursor[0]
#define xgeCURSOR_HAND      xgecursor[1]
#define xgeCURSOR_PENCIL    xgecursor[2]
#define xgeCURSOR_FLEUR     xgecursor[3]
#define xgeCURSOR_ARROW     xgecursor[4]
#define xgeCURSOR_WATCH     xgecursor[5]
#define xgeCURSOR_CIRCLE    xgecursor[6]
#define xgeCURSOR_DEFAULT   xgecursor[7]
#define xgeCURSOR_INVISIBLE xgecursor[8]

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef XGERGB_H
#include "xgergb.h"
#endif

#define xgec_MENU_BACKGROUND       xgec_Grey5
#define xgec_INFOMSG_BACKGROUND    xgec_Grey4
#define xgec_ERRORMSG_BACKGROUND   xgec_Red
#define xgec_WARNINGMSG_BACKGROUND xgec_DarkMagenta

/* ///////////////////////////////////////////////////////////////////////// */
/* wrappers around XWindow drawing and some other procedures */
#define xgeSetBackground(colour) \
  XSetBackground(xgedisplay,xgegc,colour)
#define xgeSetForeground(colour) \
  XSetForeground(xgedisplay,xgegc,colour)
#define xgeSetLineAttributes(width,line_style,cap_style,join_style) \
  XSetLineAttributes(xgedisplay,xgegc,width,line_style,cap_style,join_style)
#define xgeSetDashes(n,dash_list,offset) \
  XSetDashes(xgedisplay,xgegc,offset,dash_list,n)

#define xgeDrawRectangle(w,h,x,y) \
  XDrawRectangle(xgedisplay,xgepixmap,xgegc,x,y,w,h)
#define xgeFillRectangle(w,h,x,y) \
  XFillRectangle(xgedisplay,xgepixmap,xgegc,x,y,w,h)
#define xgeDrawString(string,x,y) \
  XDrawString(xgedisplay,xgepixmap,xgegc,x,y,string,strlen(string))
#define xgeDrawLine(x0,y0,x1,y1) \
  XDrawLine(xgedisplay,xgepixmap,xgegc,x0,y0,x1,y1)
#define xgeDrawLines(n,p) \
  XDrawLines(xgedisplay,xgepixmap,xgegc,p,n,CoordModeOrigin)
#define xgeDrawArc(w,h,x,y,a0,a1) \
  XDrawArc(xgedisplay,xgepixmap,xgegc,x,y,w,h,a0,a1)
#define xgeFillArc(w,h,x,y,a0,a1) \
  XFillArc(xgedisplay,xgepixmap,xgegc,x,y,w,h,a0,a1)
#define xgeDrawPoint(x,y) \
  XDrawPoint(xgedisplay,xgepixmap,xgegc,x,y)
#define xgeDrawPoints(n,p) \
  XDrawPoints(xgedisplay,xgepixmap,xgegc,p,n,CoordModeOrigin)
#define xgeFillPolygon(shape,n,p) \
  XFillPolygon(xgedisplay,xgepixmap,xgegc,p,n,shape,CoordModeOrigin)

#define xgeCopyRectOnScreen(w,h,x,y) \
  XCopyArea(xgedisplay,xgepixmap,xgewindow,xgegc,x,y,w,h,x,y)
#define xgeRaiseWindow() \
  XRaiseWindow(xgedisplay,xgewindow)
#define xgeResizeWindow(w,h) \
  XResizeWindow(xgedisplay,xgewindow,w,h)
#define xgeMoveWindow(x,y) \
  XMoveWindow(xgedisplay,xgewindow,x,y)

/* ///////////////////////////////////////////////////////////////////////// */
typedef struct xge_widget {
          int   id;
          short w, h, x, y;
          void *data0, *data1, *data2;
          short xofs, yofs;
          char  rpos;
          char  window_num;
          short state;
          boolean (*msgproc) ( struct xge_widget*, int, int, short, short );
          void (*redraw) ( struct xge_widget*, boolean );
          struct xge_widget *next, *prev, *up;
        } xge_widget;

typedef boolean (*xge_message_proc) ( struct xge_widget *er,
                                      int msg, int key, short x, short y );
typedef void (*xge_redraw_proc) ( struct xge_widget *er, boolean onscreen );

/* ///////////////////////////////////////////////////////////////////////// */
extern int   xge_winnum;
extern unsigned int xge_mouse_buttons;
extern int   xge_mouse_x, xge_mouse_y;
extern short xge_xx, xge_yy;

extern Display     *xgedisplay;
extern XVisualInfo *xgevisualinfo;
extern Colormap    xgecolormap;
extern int         xgescreen;
extern Window      xgeroot;
extern Window      xgewindow;
extern Pixmap      xgepixmap;
extern GC          xgegc;
extern Visual      *xgevisual;
extern XSizeHints  xgehints;
extern Cursor      xgecursor[];
extern KeySym      xgekeysym;
extern XEvent      xgeevent;

typedef struct {
         unsigned short r_bits,  g_bits,  b_bits;
         char           nr_bits, ng_bits, nb_bits;
         unsigned char  r_shift, g_shift, b_shift;
         unsigned char  r_mask, g_mask, b_mask;
         float          ar, ag, ab;
       } xge_rgbmap_bits;
extern xge_rgbmap_bits xge_rgbmap;

extern float      xge_aspect;

extern unsigned int xge_current_width, xge_current_height;  /* window dimensions */

extern char       *xge_p_name;
extern char       xge_done;
extern short      xge_prevx, xge_prevy;

extern xgecolour_int xge_foreground, xge_background;
extern int           xge_nplanes;
extern xgecolour_int *xge_palette;
extern const char    *xge_colour_name[];
extern xge_widget    *xge_null_widget;

/* ///////////////////////////////////////////////////////////////////////// */
/* widgets */
xge_widget *xge_NewWidget (
         char window_num, xge_widget *prev, int id,
         short w, short h, short x, short y,
         void *data0, void *data1,
         boolean (*msgproc)(xge_widget*, int, int, short, short),
         void (*redraw)(xge_widget*, boolean) );
void xge_SetWidgetPositioning ( xge_widget *edr,
                                char rpos, short xofs, short yofs );

void xge_DrawEmpty ( xge_widget *er, boolean onscreen );
boolean xge_EmptyMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewEmptyWidget ( char window_num, xge_widget *prev, int id,
                                 short w, short h, short x, short y );

void xge_DrawMenu ( xge_widget *er, boolean onscreen );
boolean xge_MenuMsg ( xge_widget *er, int msg, int key, short x, short y );
boolean xge_PopupMenuMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewMenu ( char window_num, xge_widget *prev, int id,
                          short w, short h, short x, short y,
                          xge_widget *widgetlist );
void xge_DrawFMenu ( xge_widget *er, boolean onscreen );
xge_widget *xge_NewFMenu ( char window_num, xge_widget *prev, int id,
                           short w, short h, short x, short y,
                           xge_widget *widgetlist );
void xge_SetMenuWidgets ( xge_widget *menu, xge_widget *widgetlist,
                          boolean redraw );

void xge_DrawSwitch ( xge_widget *er, boolean onscreen );
boolean xge_SwitchMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewSwitch ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            char *title, boolean *sw );

void xge_DrawButton ( xge_widget *er, boolean onscreen );
boolean xge_ButtonMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewButton ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            char *title );

void xge_DrawSlidebarf ( xge_widget *er, boolean onscreen );
boolean xge_SlidebarfMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewSlidebarf ( char window_num, xge_widget *prev, int id,
                               short w, short h, short x, short y,
                               float *data );

void xge_DrawSlidebarfRGB ( xge_widget *er, boolean onscreen );
xge_widget *xge_NewSlidebarfRGB ( char window_num, xge_widget *prev, int id,
                                  short w, short h, short x, short y,
                                  float *data );

void xge_DrawVSlidebarf ( xge_widget *er, boolean onscreen );
boolean xge_VSlidebarfMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewVSlidebarf ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                float *data );

void xge_DrawSlidebar2f ( xge_widget *er, boolean onscreen );
boolean xge_Slidebar2fMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewSlidebar2f ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                float *data );

void xge_DrawVSlidebar2f ( xge_widget *er, boolean onscreen );
boolean xge_VSlidebar2fMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewVSlidebar2f ( char window_num, xge_widget *prev, int id,
                                 short w, short h, short x, short y,
                                 float *data );

float xge_LinSlidebarValuef ( float xmin, float xmax, float t );
float xge_LogSlidebarValuef ( float xmin, float xmax, float t );

void xge_DrawSlidebard ( xge_widget *er, boolean onscreen );
boolean xge_SlidebardMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewSlidebard ( char window_num, xge_widget *prev, int id,
                               short w, short h, short x, short y,
                               double *data );

void xge_DrawSlidebardRGB ( xge_widget *er, boolean onscreen );
xge_widget *xge_NewSlidebardRGB ( char window_num, xge_widget *prev, int id,
                                  short w, short h, short x, short y,
                                  double *data );

void xge_DrawVSlidebard ( xge_widget *er, boolean onscreen );
boolean xge_VSlidebardMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewVSlidebard ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                double *data );

void xge_DrawSlidebar2d ( xge_widget *er, boolean onscreen );
boolean xge_Slidebar2dMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewSlidebar2d ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                double *data );

void xge_DrawVSlidebar2d ( xge_widget *er, boolean onscreen );
boolean xge_VSlidebar2dMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewVSlidebar2d ( char window_num, xge_widget *prev, int id,
                                 short w, short h, short x, short y,
                                 double *data );

double xge_LinSlidebarValued ( double xmin, double xmax, double t );
double xge_LogSlidebarValued ( double xmin, double xmax, double t );

void xge_DrawDialf ( xge_widget *er, boolean onscreen );
boolean xge_DialfMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewDialf ( char window_num, xge_widget *prev, int id,
                           short w, short h, short x, short y,
                           char *title, float *data );

void xge_DrawDiald ( xge_widget *er, boolean onscreen );
boolean xge_DialdMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewDiald ( char window_num, xge_widget *prev, int id,
                           short w, short h, short x, short y,
                           char *title, double *data );

typedef struct xge_quatrotballf {
    xge_widget  *er;
    quaternionf *q;
    trans3f     *tr;
    vector3f    axis;
    short       xc, yc, r1, r2, R;
    boolean     axis_ok, insert;
  } xge_quatrotballf;

void xge_DrawQuatRotBallf ( xge_widget *er, boolean onscreen );
boolean xge_QuatRotBallfMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewQuatRotBallf ( char window_num, xge_widget *prev, int id,
                                  short w, short h, short x, short y, short R,
                                  xge_quatrotballf *qball, quaternionf *q, trans3f *tr,
                                  char *title );

typedef struct xge_quatrotballd {
    xge_widget  *er;
    quaterniond *q;
    trans3d     *tr;
    vector3d    axis;
    short       xc, yc, r1, r2, R;
    boolean     axis_ok, insert;
  } xge_quatrotballd;

void xge_DrawQuatRotBalld ( xge_widget *er, boolean onscreen );
boolean xge_QuatRotBalldMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewQuatRotBalld ( char window_num, xge_widget *prev, int id,
                                  short w, short h, short x, short y, short R,
                                  xge_quatrotballd *qball, quaterniond *q, trans3d *tr,
                                  char *title );

void xge_DrawText ( xge_widget *er, boolean onscreen );
boolean xge_TextWidgetMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewTextWidget ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                char *text );

void xge_DrawRGBSamplef ( xge_widget *er, boolean onscreen );
xge_widget *xge_NewRGBSamplef ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                float *data );

void xge_DrawRGBSampled ( xge_widget *er, boolean onscreen );
xge_widget *xge_NewRGBSampled ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                double *data );

typedef struct xge_string_ed {
    xge_widget *er;
    short maxlength,  /* maximal string length */
          chdisp,     /* number of characters displayed */
          start,      /* first character displayed */
          pos;        /* text cursor position */
  } xge_string_ed;

void xge_DrawStringEd ( xge_widget *er, boolean onscreen );
boolean xge_StringEdMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewStringEd ( char window_num, xge_widget *prev, int id,
                              short w, short h, short x, short y,
                              short maxlength, char *text, xge_string_ed *ed );

typedef struct xge_int_widget {
    xge_widget *er;
    int        minvalue, maxvalue;
    char       *title;
  } xge_int_widget;

void xge_DrawIntWidget ( xge_widget *er, boolean onscreen );
boolean xge_IntWidgetMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewIntWidget ( char window_num, xge_widget *prev, int id,
                               short w, short h, short x, short y,
                               int minvalue, int maxvalue,
                               xge_int_widget *iw, char *title, int *valptr );

#define xge_LISTDIST 16

typedef struct xge_listbox {
    xge_widget   *er;
    char          dlistnpos;   /* number of positions displayed */
    char          maxitl;      /* maximal item length, in characters */
    short         nitems;      /* current number of list elements */
    short         fditem;      /* first displayed item */
    short         currentitem; /* current item */
    int           *itemind;    /* indexes to the item strings */
    char          *itemstr;    /* item strings */
    xgecolour_int bk0, bk1;    /* background colours for the items */
  } xge_listbox;

void xge_DrawListBox ( xge_widget *er, boolean onscreen );
boolean xge_ListBoxMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewListBox ( char window_num, xge_widget *prev, int id,
                             short w, short h, short x, short y,
                             xge_listbox *listbox );
void xge_ClearListBox ( xge_listbox *lbox );
void xge_ShortenString ( const char *s, char *buf, int maxlen );
boolean xge_GetCurrentListBoxString ( xge_listbox *lbox, char *string );
int xge_MoveInListBox ( xge_listbox *lbox, short amount );

boolean xge_SetupFileList ( xge_listbox *lbox, const char *dir,
                            const char *filter, boolean hidden );
boolean xge_SetupDirList ( xge_listbox *lbox, const char *dir,
                           const char *filter, boolean hidden,
                           const char *prevdir );
boolean xge_FilterMatches ( const char *name, const char *filter );


#define xge_2DWIN_MIN_ZOOM             0.01
#define xge_2DWIN_MAX_ZOOM           100.0

#define xge_2DWIN_NO_TOOL                0
#define xge_2DWIN_MOVING_TOOL            1
#define xge_2DWIN_SCALING_TOOL           2
#define xge_2DWIN_ROTATING_TOOL          3
#define xge_2DWIN_SHEAR_TOOL             4
#define xge_2DWIN_SELECTING_TOOL         5
#define xge_2DWIN_PANNING_TOOL           6
#define xge_2DWIN_SPECIAL_SELECTING_TOOL 7
#define xge_2DWIN_SPECIAL_TRANS_TOOL     8

typedef struct xge_2Dwinf {
    xge_widget  *er;
    CameraRecf  CPos;
    Box2f       DefBBox, RefBBox;
    boolean     panning, selecting_mode, special_selecting_mode;
    boolean     display_coord, inside;
    boolean     moving_tool, scaling_tool, rotating_tool, shear_tool,
                special_trans_tool;
    char        current_tool;
    char        tool_mode;
    int         current_tab, current_point;
    short       xx, yy;
    float       zoom;
    Box2s       selection_rect;
    point2f     saved_centre;
    point2f     scaling_centre;
    vector2f    scaling_factors;
    short       scaling_size;
    point2f     rotating_centre;
    short       rotating_radius;
    vector2f    trans_params;
    point2f     shear_centre;
    vector2f    shear_axis[2];
    float       shear_radius;
    trans2f     gwtrans;
  } xge_2Dwinf;

boolean xge_2DwinfMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_New2Dwinf ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            xge_2Dwinf *_2Dwin,
                            void (*redraw)(xge_widget*, boolean) );

void xge_2DwinfSetDefBBox ( xge_2Dwinf *_2Dwin,
                            float x0, float x1, float y0, float y1 );
void xge_2DwinfSetupProjection ( xge_2Dwinf *_2Dwin );
void xge_2DwinfPan ( xge_widget *er, short x, short y );
void xge_2DwinfZoom ( xge_widget *er, short y );
void xge_2DwinfInitProjection ( xge_2Dwinf *_2Dwin,
                                float x0, float x1, float y0, float y1 );
void xge_2DwinfResetGeomWidgets ( xge_2Dwinf *_2Dwin );
void xge_2DwinfResetGeomWidgetPos ( xge_2Dwinf *_2Dwin );
void xge_2DwinfEnableGeomWidget ( xge_2Dwinf *_2Dwin, char tool );
void xge_2DwinfDrawGeomWidgets ( xge_widget *er );
char xge_2DwinfIsItAGeomWidget ( xge_2Dwinf *_2Dwin, int key, short x, short y );
void xge_2DwinfMoveGeomWidget ( xge_2Dwinf *_2Dwin, short x, short y );
boolean xge_2DwinfApplyGeomWidget ( xge_2Dwinf *_2Dwin, short x, short y,
                                    boolean alt );
void xge_2DwinfExitWidgetMode ( xge_2Dwinf *_2Dwin );
void xge_2DwinfResetGeomWidget ( xge_2Dwinf *_2Dwin );
void xge_2DwinfDrawCursorPos ( xge_2Dwinf *_2Dwin, short x, short y );


typedef struct xge_2Dwind {
    xge_widget  *er;
    CameraRecd  CPos;
    Box2d       DefBBox, RefBBox;
    boolean     panning, selecting_mode, special_selecting_mode;
    boolean     display_coord, inside;
    boolean     moving_tool, scaling_tool, rotating_tool, shear_tool,
                special_trans_tool;
    char        current_tool;
    char        tool_mode;
    int         current_tab, current_point;
    short       xx, yy;
    double      zoom;
    Box2s       selection_rect;
    point2d     saved_centre;
    point2d     scaling_centre;
    vector2d    scaling_factors;
    short       scaling_size;
    point2d     rotating_centre;
    short       rotating_radius;
    vector2d    trans_params;
    point2d     shear_centre;
    vector2d    shear_axis[2];
    double      shear_radius;
    trans2d     gwtrans;
  } xge_2Dwind;

boolean xge_2DwindMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_New2Dwind ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            xge_2Dwind *_2Dwin,
                            void (*redraw)(xge_widget*, boolean) );

void xge_2DwindSetDefBBox ( xge_2Dwind *_2Dwin,
                            double x0, double x1, double y0, double y1 );
void xge_2DwindSetupProjection ( xge_2Dwind *_2Dwin );
void xge_2DwindPan ( xge_widget *er, short x, short y );
void xge_2DwindZoom ( xge_widget *er, short y );
void xge_2DwindInitProjection ( xge_2Dwind *_2Dwin,
                                double x0, double x1, double y0, double y1 );
void xge_2DwindResetGeomWidgets ( xge_2Dwind *_2Dwin );
void xge_2DwindResetGeomWidgetPos ( xge_2Dwind *_2Dwin );
void xge_2DwindEnableGeomWidget ( xge_2Dwind *_2Dwin, char tool );
void xge_2DwindDrawGeomWidgets ( xge_widget *er );
char xge_2DwindIsItAGeomWidget ( xge_2Dwind *_2Dwin, int key, short x, short y );
void xge_2DwindMoveGeomWidget ( xge_2Dwind *_2Dwin, short x, short y );
boolean xge_2DwindApplyGeomWidget ( xge_2Dwind *_2Dwin, short x, short y,
                                    boolean alt );
void xge_2DwindExitWidgetMode ( xge_2Dwind *_2Dwin );
void xge_2DwindResetGeomWidget ( xge_2Dwind *_2Dwin );
void xge_2DwindDrawCursorPos ( xge_2Dwind *_2Dwin, short x, short y );


typedef struct xge_fourww {
    xge_widget *er;
    xge_widget *win[4];
    float      xsfr, ysfr;
    short      splitx, splity;
    boolean    resized;
    char       zoomwin;
  } xge_fourww;

boolean xge_CompSizeFourWW ( xge_widget *er, char cs );
void xge_DrawFourWW ( xge_widget *er, boolean onscreen );
boolean xge_FourWWMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewFourWW ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            xge_widget *ww, xge_fourww *fwwdata );


#define xge_3DWIN_MIN_PARZOOM   0.01
#define xge_3DWIN_MAX_PARZOOM 100.0
#define xge_3DWIN_MIN_ZOOM      0.1
#define xge_3DWIN_MAX_ZOOM   1000.0

#define xge_3DWIN_NO_TOOL                xge_2DWIN_NO_TOOL
#define xge_3DWIN_MOVING_TOOL            xge_2DWIN_MOVING_TOOL
#define xge_3DWIN_SCALING_TOOL           xge_2DWIN_SCALING_TOOL
#define xge_3DWIN_ROTATING_TOOL          xge_2DWIN_ROTATING_TOOL
#define xge_3DWIN_SHEAR_TOOL             xge_2DWIN_SHEAR_TOOL
#define xge_3DWIN_SELECTING_TOOL         xge_2DWIN_SELECTING_TOOL
#define xge_3DWIN_PANNING_TOOL           xge_2DWIN_PANNING_TOOL
#define xge_3DWIN_SPECIAL_SELECTING_TOOL xge_2DWIN_SPECIAL_SELECTING_TOOL
#define xge_3DWIN_SPECIAL_TRANS_TOOL     xge_2DWIN_SPECIAL_TRANS_TOOL

typedef struct xge_3Dwinf {
    xge_fourww  fww;         /* this must be the first component */
    xge_widget  *cwin[4];
    CameraRecf  CPos[5];
    Box3f       DefBBox, RefBBox, WinBBox, PerspBBox;
    boolean     panning, selecting_mode, special_selecting_mode;
    boolean     display_coord;
    signed char CoordWin;
    boolean     moving_tool, scaling_tool, rotating_tool, shear_tool,
                special_trans_tool;
    char        current_tool;
    char        tool_mode;
    int         current_tab, current_point;
    short       xx, yy;
    float       perspzoom;
    Box2s       selection_rect;
    point3f     saved_centre;
    point3f     scaling_centre;
    vector3f    scaling_factors;
    short       scaling_size;
    point3f     rotating_centre;
    short       rotating_radius;
    vector3f    trans_params;
    point3f     shear_centre;
    vector3f    shear_axis[3];
    float       shear_radius;
    trans3f     gwtrans;
  } xge_3Dwinf;

xge_widget *xge_New3Dwinf ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            xge_3Dwinf *_3Dwin,
                            void (*pararedraw)(xge_widget*, boolean),
                            void (*perspredraw)(xge_widget*, boolean) );

void xge_3DwinfSetDefBBox ( xge_3Dwinf *_3Dwin, float x0, float x1,
                            float y0, float y1, float z0, float z1 );
void xge_3DwinfSetupParProj ( xge_3Dwinf *_3Dwin, Box3f *bbox );
void xge_3DwinfSetupPerspProj ( xge_3Dwinf *_3Dwin, boolean resetpos );
void xge_3DwinfUpdatePerspProj ( xge_3Dwinf *_3Dwin );
void xge_3DwinfPanParWindows ( xge_widget *er, short x, short y );
void xge_3DwinfZoomParWindows ( xge_widget *er, short y );
void xge_3DwinfPanPerspWindow ( xge_widget *er, short x, short y );
void xge_3DwinfInitProjections ( xge_3Dwinf *_3Dwin,
                   float x0, float x1, float y0, float y1, float z0, float z1 );
void xge_3DwinfResetGeomWidgets ( xge_3Dwinf *_3Dwin );
void xge_3DwinfResetGeomWidgetPos ( xge_3Dwinf *_3Dwin );
void xge_3DwinfEnableGeomWidget ( xge_3Dwinf *_3Dwin, char tool );
void xge_3DwinfDrawCursorPos ( xge_3Dwinf *_3Dwin,
                               int id, short x, short y );
void xge_3DwinfDrawSelectionRect ( xge_widget *er );
void xge_3DwinfDrawGeomWidgets ( xge_widget *er );
char xge_3DwinfIsItAGeomWidget ( xge_3Dwinf *_3Dwin, int id,
                                 int key, short x, short y );
void xge_3DwinfMoveGeomWidget ( xge_3Dwinf *_3Dwin, int id, short x, short y );
boolean xge_3DwinfApplyGeomWidget ( xge_3Dwinf *_3Dwin, int id, short x, short y,
                                    boolean alt );
void xge_3DwinfExitWidgetMode ( xge_3Dwinf *_3Dwin );
void xge_3DwinfResetGeomWidget ( xge_3Dwinf *_3Dwin );
void xge_3DwinfSavePerspCamera ( xge_3Dwinf *_3Dwin );
void xge_3DwinfSwapPerspCameras ( xge_3Dwinf *_3Dwin );


typedef struct xge_3Dwind {
    xge_fourww  fww;         /* this must be the first component */
    xge_widget  *cwin[4];
    CameraRecd  CPos[5];
    Box3d       DefBBox, RefBBox, WinBBox, PerspBBox;
    boolean     panning, selecting_mode, special_selecting_mode;
    boolean     display_coord;
    signed char CoordWin;
    boolean     moving_tool, scaling_tool, rotating_tool, shear_tool,
                special_trans_tool;
    char        current_tool;
    char        tool_mode;
    int         current_tab, current_point;
    short       xx, yy;
    double      perspzoom;
    Box2s       selection_rect;
    point3d     saved_centre;
    point3d     scaling_centre;
    vector3d    scaling_factors;
    short       scaling_size;
    point3d     rotating_centre;
    short       rotating_radius;
    vector3d    trans_params;
    point3d     shear_centre;
    vector3d    shear_axis[3];
    double      shear_radius;
    trans3d     gwtrans;
  } xge_3Dwind;

xge_widget *xge_New3Dwind ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            xge_3Dwind *_3Dwin,
                            void (*pararedraw)(xge_widget*, boolean),
                            void (*perspredraw)(xge_widget*, boolean) );

void xge_3DwindSetDefBBox ( xge_3Dwind *_3Dwin, double x0, double x1,
                            double y0, double y1, double z0, double z1 );
void xge_3DwindSetupParProj ( xge_3Dwind *_3Dwin, Box3d *bbox );
void xge_3DwindSetupPerspProj ( xge_3Dwind *_3Dwin, boolean resetpos );
void xge_3DwindUpdatePerspProj ( xge_3Dwind *_3Dwin );
void xge_3DwindPanParWindows ( xge_widget *er, short x, short y );
void xge_3DwindZoomParWindows ( xge_widget *er, short y );
void xge_3DwindPanPerspWindow ( xge_widget *er, short x, short y );
void xge_3DwindInitProjections ( xge_3Dwind *_3Dwin,
          double x0, double x1, double y0, double y1, double z0, double z1 );
void xge_3DwindResetGeomWidgets ( xge_3Dwind *_3Dwin );
void xge_3DwindResetGeomWidgetPos ( xge_3Dwind *_3Dwin );
void xge_3DwindEnableGeomWidget ( xge_3Dwind *_3Dwin, char tool );
void xge_3DwindDrawCursorPos ( xge_3Dwind *_3Dwin,
                               int id, short x, short y );
void xge_3DwindDrawSelectionRect ( xge_widget *er );
void xge_3DwindDrawGeomWidgets ( xge_widget *er );
char xge_3DwindIsItAGeomWidget ( xge_3Dwind *_3Dwin, int id,
                                 int key, short x, short y );
void xge_3DwindMoveGeomWidget ( xge_3Dwind *_3Dwin, int id, short x, short y );
boolean xge_3DwindApplyGeomWidget ( xge_3Dwind *_3Dwin, int id, short x, short y,
                                    boolean alt );
void xge_3DwindExitWidgetMode ( xge_3Dwind *_3Dwin );
void xge_3DwindResetGeomWidget ( xge_3Dwind *_3Dwin );
void xge_3DwindSavePerspCamera ( xge_3Dwind *_3Dwin );
void xge_3DwindSwapPerspCameras ( xge_3Dwind *_3Dwin );


#define xge_KNOTWIN_MIN_SCALE   0.01
#define xge_KNOTWIN_MAX_SCALE  100.0
#define xge_KNOT_EPS          1.0e-4

typedef struct {
    xge_widget    *er;
    boolean       panning, display_coord, moving_many, locked;
    boolean       closed, altkn, switchkn;
    float         akm, bkm, umin, umax, knotscf, knotshf;
    int           clcK;
    float         clcT;
    unsigned char current_mult;
    short         xx;
    int           maxknots, lastknot, degree;
    float         *knots;
    int           maxaltknots, lastaltknot, altdegree;
    float         *altknots;
    float         newknot;
    int           current_knot;
  } xge_KnotWinf;

void xge_DrawKnotWinf ( xge_widget *er, boolean onscreen );
boolean xge_KnotWinfMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewKnotWinf ( char window_num, xge_widget *prev, int id,
                              short w, short h, short x, short y,
                              xge_KnotWinf *knw, int maxknots, float *knots );

void xge_KnotWinfDrawCursorPos ( xge_KnotWinf *knw );
void xge_KnotWinfDrawAxis ( xge_KnotWinf *knw );
void xge_KnotWinfDrawKnots ( xge_KnotWinf *knw );
void xge_KnotWinfInitMapping ( xge_KnotWinf *knw, float umin, float umax );
void xge_KnotWinfZoom ( xge_KnotWinf *knw, float scf );
void xge_KnotWinfPan ( xge_KnotWinf *knw, int dxi );
void xge_KnotWinfFindMapping ( xge_KnotWinf *knw );
void xge_KnotWinfResetMapping ( xge_KnotWinf *knw );
short xge_KnotWinfMapKnot ( xge_KnotWinf *knw, float u );
float xge_KnotWinfUnmapKnot ( xge_KnotWinf *knw, short xi );
boolean xge_KnotWinfFindNearestKnot ( xge_KnotWinf *knw, int x, int y );
boolean xge_KnotWinfSetKnot ( xge_KnotWinf *knw, short x );
boolean xge_KnotWinfInsertKnot ( xge_KnotWinf *knw, short x );
boolean xge_KnotWinfRemoveKnot ( xge_KnotWinf *knw );
void xge_KnotWinfSetAltKnots ( xge_KnotWinf *knw,
               int altmaxkn, int lastaltkn, int altdeg, float *altknots );
void xge_KnotWinfSwitchAltKnots ( xge_KnotWinf *knw );


typedef struct {
    xge_widget    *er;
    boolean       panning, display_coord, moving_many, locked;
    boolean       closed, altkn, switchkn;
    double        akm, bkm, umin, umax, knotscf, knotshf;
    int           clcK;
    double        clcT;
    unsigned char current_mult;
    short         xx;
    int           maxknots, lastknot, degree;
    double        *knots;
    int           maxaltknots, lastaltknot, altdegree;
    double        *altknots;
    double        newknot;
    int           current_knot;
  } xge_KnotWind;

void xge_DrawKnotWind ( xge_widget *er, boolean onscreen );
boolean xge_KnotWindMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewKnotWind ( char window_num, xge_widget *prev, int id,
                              short w, short h, short x, short y,
                              xge_KnotWind *knw, int maxknots, double *knots );

void xge_KnotWindDrawCursorPos ( xge_KnotWind *knw );
void xge_KnotWindDrawAxis ( xge_KnotWind *knw );
void xge_KnotWindDrawKnots ( xge_KnotWind *knw );
void xge_KnotWindInitMapping ( xge_KnotWind *knw, double umin, double umax );
void xge_KnotWindZoom ( xge_KnotWind *knw, double scf );
void xge_KnotWindPan ( xge_KnotWind *knw, int dxi );
void xge_KnotWindFindMapping ( xge_KnotWind *knw );
void xge_KnotWindResetMapping ( xge_KnotWind *knw );
short xge_KnotWindMapKnot ( xge_KnotWind *knw, double u );
double xge_KnotWindUnmapKnot ( xge_KnotWind *knw, short xi );
boolean xge_KnotWindFindNearestKnot ( xge_KnotWind *knw, int x, int y );
boolean xge_KnotWindSetKnot ( xge_KnotWind *knw, short x );
boolean xge_KnotWindInsertKnot ( xge_KnotWind *knw, short x );
boolean xge_KnotWindRemoveKnot ( xge_KnotWind *knw );
void xge_KnotWindSetAltKnots ( xge_KnotWind *knw,
               int altmaxkn, int lastaltkn, int altdeg, double *altknots );
void xge_KnotWindSwitchAltKnots ( xge_KnotWind *knw );


typedef struct {
    xge_widget  *er;
    CameraRecf  CPos;
    Box2f       DefBBox, RefBBox;
    point3f     centre;
    boolean     panning, selecting_mode, moving_many, locked_u, locked_v;
    boolean     display_coord, inside;
    unsigned char current_mult;
    int         current_item;
    short       knot_margin;
    short       xx, yy;
    Box2s       selection_rect;
    boolean     closed_u, closed_v;
    int         clcKu, clcKv;
    float       clcTu, clcTv;
    int         maxknots_u, lastknot_u, degree_u;
    float       *knots_u;
    int         maxknots_v, lastknot_v, degree_v;
    float       *knots_v;
    float       newknot;
    boolean     altknu, switchknu, altknv, switchknv;
    int         altmaxknots_u, altlastknot_u, altdeg_u;
    float       *altknots_u;
    int         altmaxknots_v, altlastknot_v, altdeg_v;
    float       *altknots_v;
  } xge_T2KnotWinf;

boolean xge_T2KnotWinfMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewT2KnotWinf ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                short knot_margin,
                                xge_T2KnotWinf *T2win,
                                void (*redraw)(xge_widget*, boolean),
                                int maxknots_u, float *knots_u,
                                int maxknots_v, float *knots_v );

void xge_T2KnotWinfDrawKnots ( xge_T2KnotWinf *T2win );
void xge_T2KnotWinfSetupMapping ( xge_T2KnotWinf *T2win );
void xge_T2KnotWinfInitMapping ( xge_T2KnotWinf *T2win,
                                 float umin, float umax, float vmin, float vmax );
void xge_T2KnotWinfZoom ( xge_T2KnotWinf *T2win, short y );
boolean xge_T2KnotWinfPan ( xge_T2KnotWinf *T2win, short x, short y );
void xge_T2KnotWinfFindMapping ( xge_T2KnotWinf *T2win );
void xge_T2KnotWinfResetMapping ( xge_T2KnotWinf *T2win );

char xge_T2KnotWinfFindDomWinRegion ( xge_T2KnotWinf *T2win, int x, int y );
char xge_T2KnotWinfFindNearestKnot ( xge_T2KnotWinf *T2win, int x, int y );
short xge_T2KnotWinfMapKnotU ( xge_T2KnotWinf *T2win, float u );
float xge_T2KnotWinfUnmapKnotU ( xge_T2KnotWinf *T2win, short xi );
short xge_T2KnotWinfMapKnotV ( xge_T2KnotWinf *T2win, float v );
float xge_T2KnotWinfUnmapKnotV ( xge_T2KnotWinf *T2win, short eta );
boolean xge_T2KnotWinfSetKnotU ( xge_T2KnotWinf *T2win, short x );
boolean xge_T2KnotWinfInsertKnotU ( xge_T2KnotWinf *T2win, short x );
boolean xge_T2KnotWinfRemoveKnotU ( xge_T2KnotWinf *T2win );
boolean xge_T2KnotWinfSetKnotV ( xge_T2KnotWinf *T2win, short y );
boolean xge_T2KnotWinfInsertKnotV ( xge_T2KnotWinf *T2win, short y );
boolean xge_T2KnotWinfRemoveKnotV ( xge_T2KnotWinf *T2win );
void xge_T2KnotWinfSelect ( xge_T2KnotWinf *T2win,
                            short x0, short x1, short y0, short y1 );
void xge_T2KnotWinfUnselect ( xge_T2KnotWinf *T2win,
                              short x0, short x1, short y0, short y1 );

void xge_T2KnotWinfSetAltKnots ( xge_T2KnotWinf *T2win,
               int altmaxknu, int lastaltknu, int altdegu, float *altknotsu,
               int altmaxknv, int lastaltknv, int altdegv, float *altknotsv );
void xge_T2KnotWinfSwitchAltKnots ( xge_T2KnotWinf *T2win,
               boolean altu, boolean altv );
void xge_T2KnotWinfDrawCursorPos ( xge_T2KnotWinf *T2win, short x, short y );


typedef struct {
    xge_widget  *er;
    CameraRecd  CPos;
    Box2d       DefBBox, RefBBox;
    point3d     centre;
    boolean     panning, selecting_mode, moving_many, locked_u, locked_v;
    boolean     display_coord, inside;
    unsigned char current_mult;
    int         current_item;
    short       knot_margin;
    short       xx, yy;
    Box2s       selection_rect;
    boolean     closed_u, closed_v;
    int         clcKu, clcKv;
    double      clcTu, clcTv;
    int         maxknots_u, lastknot_u, degree_u;
    double      *knots_u;
    int         maxknots_v, lastknot_v, degree_v;
    double      *knots_v;
    double      newknot;
    boolean     altknu, switchknu, altknv, switchknv;
    int         altmaxknots_u, altlastknot_u, altdeg_u;
    double      *altknots_u;
    int         altmaxknots_v, altlastknot_v, altdeg_v;
    double      *altknots_v;
  } xge_T2KnotWind;

boolean xge_T2KnotWindMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewT2KnotWind ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                short knot_margin,
                                xge_T2KnotWind *T2win,
                                void (*redraw)(xge_widget*, boolean),
                                int maxknots_u, double *knots_u,
                                int maxknots_v, double *knots_v );

void xge_T2KnotWindDrawKnots ( xge_T2KnotWind *T2win );
void xge_T2KnotWindSetupMapping ( xge_T2KnotWind *T2win );
void xge_T2KnotWindInitMapping ( xge_T2KnotWind *T2win,
                                 double umin, double umax, double vmin, double vmax );
void xge_T2KnotWindZoom ( xge_T2KnotWind *T2win, short y );
boolean xge_T2KnotWindPan ( xge_T2KnotWind *T2win, short x, short y );
void xge_T2KnotWindFindMapping ( xge_T2KnotWind *T2win );
void xge_T2KnotWindResetMapping ( xge_T2KnotWind *T2win );

char xge_T2KnotWindFindDomWinRegion ( xge_T2KnotWind *T2win, int x, int y );
char xge_T2KnotWindFindNearestKnot ( xge_T2KnotWind *T2win, int x, int y );
short xge_T2KnotWindMapKnotU ( xge_T2KnotWind *T2win, double u );
double xge_T2KnotWindUnmapKnotU ( xge_T2KnotWind *T2win, short xi );
short xge_T2KnotWindMapKnotV ( xge_T2KnotWind *T2win, double v );
double xge_T2KnotWindUnmapKnotV ( xge_T2KnotWind *T2win, short eta );
boolean xge_T2KnotWindSetKnotU ( xge_T2KnotWind *T2win, short x );
boolean xge_T2KnotWindInsertKnotU ( xge_T2KnotWind *T2win, short x );
boolean xge_T2KnotWindRemoveKnotU ( xge_T2KnotWind *T2win );
boolean xge_T2KnotWindSetKnotV ( xge_T2KnotWind *T2win, short y );
boolean xge_T2KnotWindInsertKnotV ( xge_T2KnotWind *T2win, short y );
boolean xge_T2KnotWindRemoveKnotV ( xge_T2KnotWind *T2win );
void xge_T2KnotWindSelect ( xge_T2KnotWind *T2win,
                            short x0, short x1, short y0, short y1 );
void xge_T2KnotWindUnselect ( xge_T2KnotWind *T2win,
                              short x0, short x1, short y0, short y1 );

void xge_T2KnotWindSetAltKnots ( xge_T2KnotWind *T2win,
               int altmaxknu, int lastaltknu, int altdegu, double *altknotsu,
               int altmaxknv, int lastaltknv, int altdegv, double *altknotsv );
void xge_T2KnotWindSwitchAltKnots ( xge_T2KnotWind *T2win,
               boolean altu, boolean altv );
void xge_T2KnotWindDrawCursorPos ( xge_T2KnotWind *T2win, short x, short y );


typedef struct {
    xge_widget *er;
    xge_widget *contents, *clipw, *xsl, *ysl;
    float      x, y;
    boolean    xslon, yslon;
  } xge_scroll_widget;

void xge_SetupScrollWidgetPos ( xge_widget *er );
void xge_DrawScrollWidget ( xge_widget *er, boolean onscreen );
boolean xge_ScrollWidgetMsg ( xge_widget *er, int msg, int key, short x, short y );
xge_widget *xge_NewScrollWidget ( char window_num, xge_widget *prev, int id,
                                  short w, short h, short x, short y,
                                  xge_scroll_widget *sw, xge_widget *contents );

/* ///////////////////////////////////////////////////////////////////////// */
void xge_AddPopup ( xge_widget *er );
void xge_RemovePopup ( boolean redraw );
void xge_RemovePopups ( boolean redraw );
boolean xge_IsPopupOn ( xge_widget *er );

/* ///////////////////////////////////////////////////////////////////////// */
/* error and info message procedures */
void _xge_DisplayErrorMessage ( char *message, xgecolour_int bk, int key );
void xge_DisplayErrorMessage ( char *message, int key );
void xge_DisplayWarningMessage ( char *message, int key );
void xge_DisplayInfoMessage ( char **msglines, int key );

/* ///////////////////////////////////////////////////////////////////////// */
void xge_OutPixels ( const xpoint *buf, int n );
boolean xge_DrawBC2f ( int n, const point2f *cp );
boolean xge_DrawBC2Rf ( int n, const point3f *cp );
boolean xge_DrawBC2d ( int n, const point2d *cp );
boolean xge_DrawBC2Rd ( int n, const point3d *cp );

int     xge_NewWindow ( char *title );
boolean xge_SetWindow ( int win );
int     xge_CurrentWindow ( void );
void    xge_SetWinEdRect ( xge_widget *edr );

int  xge_NewCursor ( int shape );
void xge_SetWindowCursor ( int win, Cursor cursor );
void xge_SetCurrentWindowCursor ( Cursor cursor );
void xge_SetOtherWindowsCursor ( Cursor cursor );

void xge_RedrawPopups ( void );
void xge_Redraw ( void );
void xge_RedrawAll ( void );

boolean xge_PointInRect ( xge_widget *edr, short x, short y );
void xge_BoundPoint ( xge_widget *er, short *x, short *y );
boolean xge_RectanglesIntersect ( short wa, short ha, short xa, short ya,
                                  short wb, short hb, short xb, short yb );

boolean xge_IntersectXRectangles ( XRectangle *r1, XRectangle *r2 );
boolean xge_SetClipping ( xge_widget *edr );
void xge_ResetClipping ( void );

void xge_RepositionWidgets ( short w, short h, short x, short y,
                             xge_widget *edr );

xgecolour_int xge_PixelColourf ( float r, float g, float b );
xgecolour_int xge_PixelColour ( byte r, byte g, byte b );
void xge_GetPixelColour ( xgecolour_int pixel, byte *r, byte *g, byte *b );
void xge_GetPixelColourf ( xgecolour_int pixel, float *r, float *g, float *b );

void xge_DrawVShadedRect ( short w, short h, short x, short y,
                           xgecolour_int c0, xgecolour_int c1, short nb );
void xge_DrawHShadedRect ( short w, short h, short x, short y,
                           xgecolour_int c0, xgecolour_int c1, short nb );

void xge_OrderSelectionRect ( Box2s *sel_rect );
void xge_DrawGeomWinBackground ( xge_widget *er );
void xge_DrawGeomWinFrame ( xge_widget *er, boolean onscreen );
void xge_DrawGeomWinSelectionRect ( xge_widget *er, Box2s *sel_rect );

void xge_GrabFocus ( xge_widget *er, boolean all );
void xge_ReleaseFocus ( xge_widget *er );
xge_widget *xge_GetFocusWidget ( char win );

void xge_GetWindowSize ( void );
void xge_MessageLoop ( void );
void xge_PostIdleCommand ( unsigned int key, short x, short y );
void xge_FlushFromAThread ( void );

void xge_dispatch_message ( unsigned int msg, unsigned int key, short x, short y );
void xge_get_message ( unsigned int *msg, unsigned int *key, short *x, short *y );

/* ///////////////////////////////////////////////////////////////////////// */
/* The following procedures are intended for pure Xlib applications.         */
/* In OpenGL (GLX) application use xgle_... procedures, whose headers are    */
/* in the file xgledit.h.                                                    */
#ifndef XGLEDIT_H
void xge_Init ( int argc, char *argv[],
                int (*callback)(xge_widget*,int,int,short,short),
                char *title );
void xge_Cleanup ( void );
#endif

#ifdef __cplusplus
}
#endif

#endif /*XGEDIT_H*/

