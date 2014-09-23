
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2010                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libpsout library of C procedures -                */
/* generating PostScript (TM) files                                      */ 

#ifndef PSOUT_H
#define PSOUT_H

#ifndef PKGEOM_H
#include "pkgeom.h"
#endif

#define PS_FILE_STACK_LENGTH 4

#ifdef __cplusplus   
extern "C" {
#endif

extern short ps_dec_digits;

boolean ps_OpenFile ( const char *filename, unsigned int dpi );
void ps_CloseFile ( void );

void ps_Set_Gray ( float gray );
void ps_Set_RGB ( float red, float green, float blue );
void ps_Set_Line_Width ( float w );

void ps_Draw_Line ( float x1, float y1, float x2, float y2 );
void ps_Set_Clip_Rect ( float w, float h, float x, float y );
void ps_Draw_Rect ( float w, float h, float x, float y );
void ps_Fill_Rect ( float w, float h, float x, float y );
void ps_Hatch_Rect ( float w, float h, float x, float y,
                     float ang, float d );

void ps_Draw_Polyline2f ( int n, const point2f *p );
void ps_Draw_Polyline2d ( int n, const point2d *p );
void ps_Draw_Polyline2Rf ( int n, const point3f *p );
void ps_Draw_Polyline2Rd ( int n, const point3d *p );
void ps_Set_Clip_Polygon2f ( int n, const point2f *p );
void ps_Set_Clip_Polygon2d ( int n, const point2d *p );
void ps_Set_Clip_Polygon2Rf ( int n, const point3f *p );
void ps_Set_Clip_Polygon2Rd ( int n, const point3d *p );
void ps_Fill_Polygon2f ( int n, const point2f *p );
void ps_Fill_Polygon2d ( int n, const point2d *p );
void ps_Fill_Polygon2Rf ( int n, const point3f *p );
void ps_Fill_Polygon2Rd ( int n, const point3d *p );

void ps_Draw_BezierCf ( const point2f *p, int n );
void ps_Draw_BezierCd ( const point2d *p, int n );
void ps_Draw_Circle ( float x, float y, float r );
void ps_Fill_Circle ( float x, float y, float r );
void ps_Mark_Circle ( float x, float y );
void ps_Draw_Arc ( float x, float y, float r, float a0, float a1 );
void ps_Write_Command ( char *command );

void ps_Init_Bitmap ( int w, int h, int x, int y, byte b );
void ps_Init_BitmapP ( int w, int h, int x, int y );
void ps_Init_BitmapRGB ( int w, int h, int x, int y );
void ps_Init_BitmapRGBP ( int w, int h, int x, int y );
void ps_Out_Line ( byte *data );
void ps_Out_LineP ( byte *data );
void ps_Out_LineRGB ( byte *data );
void ps_Out_LineRGBP ( byte *data );

void ps_Newpath ( void );
void ps_MoveTo ( float x, float y );
void ps_LineTo ( float x, float y );
void ps_ShCone ( float x, float y, float x1, float y1, float x2, float y2 );

void ps_GSave ( void );
void ps_GRestore ( void );

void ps_DenseScreen ( void );
void ps_WriteBBox ( float x1, float y1, float x2, float y2 );
void ps_GetSize ( float *x1, float *y1, float *x2, float *y2 );

void ps_BeginDict ( int n );
void ps_EndDict ( void );

#define tickl  10.0
#define tickw   2.0
#define tickd   6.0
#define dotr   12.0
#define arrowl 71.0
#define arroww 12.5

void psl_SetLine ( float x1, float y1, float x2, float y2, float t1,
                   float t2 );
void psl_GetPointf ( float t, float *x, float *y );
void psl_GetPointd ( double t, double *x, double *y );
float psl_GetDParam ( float dl );
void psl_GoAlong ( float s, float *x, float *y );
void psl_GoPerp ( float s, float *x, float *y );

void psl_Tick ( float t );
void psl_BTick ( float t );
void psl_HTick ( float t, boolean left );
void psl_Dot ( float t );
void psl_HDot ( float t );
void psl_TrMark ( float x, float y );
void psl_BlackTrMark ( float x, float y );
void psl_HighTrMark ( float x, float y );
void psl_BlackHighTrMark ( float x, float y );
void psl_LTrMark ( float t );
void psl_BlackLTrMark ( float t );
void psl_HighLTrMark ( float t );
void psl_BlackHighLTrMark ( float t );
void psl_Arrow ( float t, boolean sgn );
void psl_BkArrow ( float t, boolean sgn );
void psl_Draw ( float ta, float tb, float w );
void psl_ADraw ( float ta, float tb, float ea, float eb, float w );
void psl_MapsTo ( float t );
void psl_BkMapsTo ( float t );
void psl_DrawEye ( float t, byte cc, float mag, float ang );

#ifdef __cplusplus
}
#endif

#endif
