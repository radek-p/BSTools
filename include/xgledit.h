
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef XGLEDIT_H
#define XGLEDIT_H

#ifndef _LIBC_LIMITS_H_
#include <limits.h>
#endif

#ifndef __gl_h_
#include <GL/gl.h>
#endif
#ifndef __glx_h__
#include <GL/glx.h>
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
#ifndef XGEDIT_H
#include "xgedit.h"
#endif
#ifndef XGLERGB_H
#include "xglergb.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define xgleCopyGLRect(w,h,x,y) \
  XCopyArea(xgedisplay,xglepixmap,xgepixmap,xgegc,0,xge_MAX_HEIGHT-h,w,h,x,y)
#define xgleClearColor3fv(rgb) \
  glClearColor(rgb[0],rgb[1],rgb[2],1.0)

extern Pixmap      xglepixmap;
extern XID         _xglepixmap;
extern void        *xglecontext;

extern GLfloat xgle_palette[XGLE_PALETTE_LENGTH][3];

void xgle_Init ( int argc, char *argv[],
                 int (*callback)(xge_widget*,int,int,short,short),
                 char *title,
                 boolean depth, boolean accum, boolean stencil );
void xgle_Cleanup ( void );

void xgle_SetIdentMapping ( xge_widget *er );

boolean xgle_SetGLCameraf ( CameraRecf *CPos );
boolean xgle_SetGLCamerad ( CameraRecd *CPos );

boolean xgle_SetGLaccCameraf ( CameraRecf *CPos,
                               float pixdx, float pixdy,
                               float eyedx, float eyedy, float focus );
boolean xgle_SetGLaccCamerad ( CameraRecd *CPos,
                               double pixdx, double pixdy,
                               double eyedx, double eyedy, double focus );

void xgle_MultMatrix3f ( trans3f *tr );
void xgle_MultMatrix3d ( trans3d *tr );

void xgleDrawPoint ( int x, int y );
void xgleDrawPoints ( int n, XPoint *p );
void xgleDrawLine ( int x0, int y0, int x1, int y1 );
void xgleDrawRectangle ( int w, int h, int x, int y );
void xgleDrawString ( char *string, int x, int y );


void xgle_DrawGeomWinBackground ( xge_widget *er, GLbitfield mask );

void xgle_2DwinfDrawCursorPos ( xge_2Dwinf *_2Dwin, short x, short y );
void xgle_2DwinfDrawAxes ( xge_2Dwinf *_2Dwin );
void xgle_2DwindDrawCursorPos ( xge_2Dwind *_2Dwin, short x, short y );
void xgle_2DwindDrawAxes ( xge_2Dwind *_2Dwin );

void xgle_3DwinfDrawCursorPos ( xge_3Dwinf *_3Dwin,
                                int id, short x, short y );
void xgle_3DwindDrawCursorPos ( xge_3Dwind *_3Dwin,
                                int id, short x, short y );

void xgle_T2KnotWinfDrawCursorPos ( xge_T2KnotWinf *T2win, short x, short y );
void xgle_T2KnotWindDrawCursorPos ( xge_T2KnotWind *T2win, short x, short y );

#ifdef __cplusplus
}
#endif

#endif /*XGLEDIT_H*/

