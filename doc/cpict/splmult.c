
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h> 
#include <string.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "psout.h"

char fn[] = "splmult.ps";

/* ///////////////////////////////////////////////////////////////////////// */
static void DrawBSPolyline ( int degree, int lastknot, const point2f *d, int mk )
{
  int i;

  ps_Draw_Polyline2f ( lastknot-degree, d );

  switch ( mk ) {
case 0:
    for ( i = 0; i < lastknot-degree; i++ )
      ps_Mark_Circle ( d[i].x, d[i].y );
    break;
case 1:
    for ( i = 0; i < lastknot-degree; i++ )
      ps_Fill_Circle ( d[i].x, d[i].y, 8.0 );
    break;
default:
    ;
  }
} /*DrawBSPolyline*/

static void DrawBSCurve ( int degree, int lastknot,
                          const float *knots, const point2f *d )
{
#define DD 200
  int i;
  point2f c[DD+1];
  float t;

  for ( i = 0; i <= DD; i++ ) {
    t = knots[degree]+(float)i/(float)DD*(knots[lastknot-degree]-knots[degree]);
    mbs_deBoorC2f ( degree, lastknot, knots, d, t, &c[i] );
  }
  ps_Draw_Polyline2f ( DD+1, c );
#undef DD
} /*DrawBSCurve*/

static void DrawGraphs ( int lastknot, const float *knots,
                         float ax, float bx, float by )
{
  int     i;
  float   t, x, y;

  ps_Set_Gray ( 0.0 );
  psl_SetLine ( bx, by, bx+ax, by, knots[0], knots[lastknot] );
  psl_Draw ( knots[0], (float)(knots[lastknot]+1.0), 2.0 );
  psl_Arrow ( (float)(knots[lastknot]+1.1), true );
  t = (float)(knots[0]-1.0e10);
  for ( i = 0; i <= lastknot; i++ ) {
    psl_Tick ( knots[i] );
    if ( knots[i] > t ) {
      psl_GetPointf ( knots[i], &x, &y );
    }
    y -= 40.0;
    ps_Fill_Circle ( x, y, 10.0 );
    t = knots[i];
  }
} /*DrawGraphs*/

static void DrawGraph1 ( int degree, int lastknot, const float *knots,
                         const float *coef,
                         float ax, float bx, float ay, float by )
{
#define DD 200
  int     i;
  float   t, x, y, f;
  point2f c[DD+1];

  ps_Set_Gray ( 0.0 );
  psl_SetLine ( bx, by, bx, by+ay, 0.0, 1.0 );
  psl_Draw ( 0.0, 1.4, 2.0 );
  psl_Tick ( 0.0 );
  psl_Tick ( 1.0 );
  psl_GetPointf ( 1.0, &x, &y );
  psl_Arrow ( 1.45, true );
  psl_SetLine ( x, y, x+ax, y, 0.0, 1.0 );
  psl_Draw ( 0.0, 1.0, 1.4 );
  psl_SetLine ( bx, by, bx+ax, by, knots[0], knots[lastknot] );
  psl_Draw ( knots[0], (float)(knots[lastknot]+1.0), 2.0 );
  psl_Arrow ( (float)(knots[lastknot]+1.1), true );
  t = (float)(knots[0]-1.0e10);
  for ( i = 0; i <= lastknot; i++ ) {
    psl_Tick ( knots[i] );
    if ( knots[i] > t ) {
      psl_GetPointf ( knots[i], &x, &y );
    }
    y -= 40.0;
    ps_Fill_Circle ( x, y, 10.0 );
    t = knots[i];
  }
  ps_Set_Gray ( 0.68 );
  for ( i = 0; i <= DD; i++ ) {
    t = knots[degree]+(float)i/(float)DD*(knots[lastknot-degree]-knots[degree]);
    psl_GetPointf ( t, &c[i].x, &c[i].y );
    mbs_deBoorC1f ( degree, lastknot, knots, coef, t, &f );
    c[i].y += ay*f;
  }
  ps_Set_Gray ( 0.68 );
  ps_Set_Line_Width ( 6.0 );
  ps_Draw_Polyline2f ( DD+1, c );
#undef DD
} /*DrawGraphs*/

/* ///////////////////////////////////////////////////////////////////////// */

#define DEGREE    3
#define LASTUKNOT 10
float u[LASTUKNOT+1] = {0.0, 1.0, 1.0, 1.0, 4.0, 4.0, 7.0, 9.0, 9.0, 9.0, 10.0};
point2f d[LASTUKNOT-DEGREE] =
  {{250.0,2100.0},{250.0,2400.0},{400.0,2500.0},{800.0,2400.0},
   {900.0,2600.0},{1200.0,2500.0},{1700.0,2100.0}};

#define DEGS      2
#define LASTVKNOT 6
float   v[LASTVKNOT+1] = {0.0, 1.0, 1.0, 4.5, 9.0, 9.0, 10.0};
float   s[LASTVKNOT-DEGS] = {1.0, 2.0, 0.0, 1.6};

int     outdegree, outlastkn;
float   outkn[100];
point2f outcp[100];

#define CX 600.0
#define CY 1950.0

static void Exp1 ()
{
  int i;

  ps_WriteBBox ( 10, 83, 304, 337 );
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();

  outlastkn = mbs_BSProdRepSizef ( DEGS, LASTVKNOT, v, DEGREE, LASTUKNOT, u );
  for ( i = 0; i < LASTUKNOT-DEGREE; i++ ) {
    d[i].x -= CX;
    d[i].y -= CY;
  }
  mbs_multiMultBSCf ( 1, DEGS, LASTVKNOT, v, 0, s,
                      2, 1, DEGREE, LASTUKNOT, u, 0, (float*)d, &outdegree,
                      &outlastkn, outkn, 0, (float*)outcp );
  for ( i = 0; i < LASTUKNOT-DEGREE; i++ ) {
    d[i].x += CX;
    d[i].y += CY;
  }
  for ( i = 0; i < outlastkn-outdegree; i++ ) {
    outcp[i].x += CX;
    outcp[i].y += CY;
  }
  ps_Set_Gray ( 0.0 );
  psl_SetLine ( CX, CY, d[0].x, d[0].y, 0.0, 1.0 );
  psl_Draw ( 0.0, 0.95, 2.0 );
  psl_Arrow ( 1.0, true );
  psl_SetLine ( CX, CY, d[LASTUKNOT-DEGREE-1].x, d[LASTUKNOT-DEGREE-1].y, 0.0, 1.0 );
  psl_Draw ( 0.0, s[LASTVKNOT-DEGS-1], 2.0 );
  psl_Arrow ( 1.0, true );
  psl_Arrow ( s[LASTVKNOT-DEGS-1], true );
  ps_Set_Line_Width ( 1.4 );
  ps_Draw_Line ( CX, (float)(CY-150.0), CX, (float)(CY+150.0) );
  ps_Draw_Line ( (float)(CX-150.0), CY, (float)(CX+150.0), CY );
  ps_Mark_Circle ( CX, CY );
  ps_Set_Line_Width ( 8.0 );
  DrawBSCurve ( DEGREE, LASTUKNOT, u, d );
  ps_Set_Gray ( 0.72 );
  ps_Set_Line_Width ( 2.0 );
  DrawBSCurve ( outdegree, outlastkn, outkn, outcp );

  ps_Set_Line_Width ( 4.0 );
  ps_Set_Gray ( 0.0 );
  DrawBSPolyline ( DEGREE, LASTUKNOT, d, 0 );

  ps_Set_Line_Width ( 2.0 );
  ps_Set_Gray ( 0.68 );
  DrawBSPolyline ( outdegree, outlastkn, outcp, 1 );

  DrawGraph1 ( DEGS, LASTVKNOT, v, s, 2200.0, 100.0, 300.0, 1400.0 );
  DrawGraphs ( LASTUKNOT, u, 2200.0, 100.0, 1150.0 );
  DrawGraphs ( outlastkn, outkn, 2200.0, 100.0, 900.0 );

  ps_CloseFile ();
  printf ( "%s\n", fn );
} /*Exp1*/

int main ()
{
  pkv_InitScratchMem ( 65536 );
  Exp1 ();
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

