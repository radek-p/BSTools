
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* the definition below  seems to be completely irrelevant, window managers */
/* ignore it */
#define xge_BORDER_WIDTH 2

/* xgedit X window descriptor */
typedef struct {
    Window     thewindow;
    Pixmap     thepixmap;
    int        win_rect_num;
    xge_widget *win_edr0, *win_edr1;
    xge_widget *popup0, *popup1;
    Cursor     cursor;
    Cursor     cursorstack[xge_FOCUS_DEPTH];
    xge_widget *focusstack[xge_FOCUS_DEPTH];
    char       fsp;  /* focus stack pointer */
    XRectangle thewinrect;
  } WinDesc;

extern WinDesc xge_windesc[xge_MAX_WINDOWS];  /* window descriptors */

extern int     xge_current_win;
extern boolean xge_notinfocus;

extern int _xge_argc;
extern char **_xge_argv;

extern xge_widget    *xge_lastwin;
extern int           xge_usr_msg_key;
extern xge_widget    *xge_errmsg_edr;
extern int           xge_errmsg_win;
extern char          *xge_errmsg_msgtext;
extern char          **xge_info_msgtext;
extern xgecolour_int xge_msgbkcolour;

extern XID xgeglxpixmap;

extern int (*xge_callback)(xge_widget*,int,int,short,short);

extern short xge_cursnum;

extern xge_widget *_xge_background_widget;
extern xge_widget *_xge_special_widget;

boolean _xge_background_msg ( xge_widget *er,
                             int msg, int key, short x, short y );
void _xge_FindAspect ( void );

void xge_RemoveErrorMessage ( void );

boolean xge_InitRectAllocation ( void );
xge_widget *xge_AllocEdRect ( void );
void xge_FreeEdRectangles ( void );

void xge_SetupEdRect ( char window_num, xge_widget *edr, int en, int id,
                       short w, short h, short x, short y,
                       boolean (*msgproc) ( xge_widget*, int, int, short, short ),
                       void (*redraw) ( xge_widget*, boolean ) );

#ifdef XGLEDIT_H
void _xgle_MakePalette ( void );
#else
void _xge_MakePalette ( void );
#endif

boolean xge_CallMsgProc ( xge_widget **last, xge_widget *er,
                          int msg, int key, short x, short y );

void _xge_QuatRotBallDrawCircles ( short xc, short yc, short r, trans3f *tr,
                                   void (*outpixels)(const xpoint *buf, int n) );

#ifdef USE_XEXT_SHAPE
/* ///////////////////////////////////////////////////////////////////////// */
/* If X Window shape extension is available, widgets may be drawn using      */
/* a special window, which may actually stick beyond the area of the regular */
/* application window. The window has no border and no title bar.            */
typedef struct {
    Window     thewindow;
    XRectangle thewinrect;
    short      xpos, ypos, xpix, ypix;
    boolean    nonempty, mapped;
  } SpecialWinDesc;

extern boolean        xge_try_ext, xge_use_specialwin, xge_specialwin_in_use;
extern SpecialWinDesc xge_specialwin[4];
extern Pixmap         xge_specialpixmap;
extern XRectangle     xge_specialpixmaprect;
extern GC             xge_specialwingc, xge_specialpixmapgc;
extern int            xge_specialevent_base, xge_specialerror_base;

extern void (*_xge_compspecialwinsizes)( xge_widget *wdg );
extern void (*_xge_maskspecialwin)(xge_widget *wdg);
extern void (*_xge_drawspecialwin)(int w, xge_widget *wdg);
extern xge_widget *_xge_specialwdg;

boolean _xge_CreateSpecialWin ( void );
boolean _xge_MapSpecialWin ( int trwin,
                             void (*compspecialwinsizes)(xge_widget *wdg),
                             void (*maskspecialwin)(xge_widget *wdg),
                             void (*drawspecialwin)(int w, xge_widget *wdg),
                             xge_widget *wdg );
void _xge_UnmapSpecialWin ( void );
void _xge_DestroySpecialWin ( void );
void _xge_RemaskSpecialWin ( void );
void _xge_SpecialWinEvent ( void );
void _xge_OutSpecialMaskPixels ( const xpoint *buf, int n );

void _xge_QuatRotCompSpecialWinSizes ( xge_widget *spqw );
void _xge_QuatRotBallDrawSpecialWin ( int w, xge_widget *wdg );
#endif
