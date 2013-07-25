
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#define xge_BORDER_WIDTH 2

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

void _xge_QuatRotBallDrawCircles ( short xc, short yc, short r, trans3f *tr );

