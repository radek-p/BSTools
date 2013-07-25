
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "psout.h"

char fn[] = "bsder123.ps";

static void DrawCPolygon ( int degree, int lastknot, const point2f *d )
{
  int i;

  ps_Draw_Polyline2f ( lastknot-degree, d );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < lastknot-degree; i++ )
    ps_Mark_Circle ( d[i].x, d[i].y );
} /*DrawCPolygon*/

static void DrawBSCurve ( int degree, int lastknot, const float *u,
                          const point2f *d )
{
#define DD 200
  int     i;
  float   t;
  point2f *c;

  c = pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  for ( i = 0; i <= DD; i++ ) {
    t = u[degree] + (float)i/(float)DD*(u[lastknot-degree]-u[degree]);
    mbs_deBoorC2f ( degree, lastknot, u, d, t, &c[i] );
  }
  ps_Draw_Polyline2f ( DD+1, c );
  pkv_FreeScratchMem ( (DD+1)*sizeof(point2f) );
#undef DD
} /*DrawBSCurve*/

static void DrawKnots ( int NNN, float *kn, char c )
{
  int i;
  float w, x, y;

  if ( c ) {
    ps_Set_Gray ( 0.0 );  w = 2.0;
  }
  else {
    ps_Set_Gray ( 0.6 );  w = 6.0;
  }
  psl_SetLine ( 1600.0, 350.0, 3000.0, 350.0, kn[1], kn[NNN-1] );
  psl_Draw ( kn[1], kn[NNN-1], w );
  for ( i = 1; i < NNN; i++ ) {
    if ( i == 1 || kn[i] > kn[i-1] ) {
      psl_Tick ( kn[i] );
      psl_GetPointf ( kn[i], &x, &y );
    }
    y -= 40;
    ps_Fill_Circle ( x, y, 10 );
  }
} /*DrawKnots*/

/* ///////////////////////////////////////////////////////////////////////// */
#define N 15
#define n  5
float u[N+1] = {0.0,0.0,0.0,0.0,0.0,0.0,
                3.5,3.5,3.5,3.5,
                7.0,7.0,7.0,7.0,7.0,7.0};

point2f d[N-n] =
  {{640.0,400.0},{320.0,400.0},{160.0,720.0},{640.0,1200.0},{1120.0,1200.0},
   {1440.0,880.0},{1120.0,720.0},{1120.0,400.0},{1440.0,400.0},{1600.0,560.0}};

point2f  p;
vector2f d1, d2, d3;

int main ( void )
{
#define DD 60
  int   i;
  float t;

  pkv_InitScratchMem ( 65536 );

  mbs_multideBoorDer3f ( n, N, u, 1, 2, 0, (float*)d,
                       3.5, (float*)&p, (float*)&d1, (float*)&d2, (float*)&d3 );
/*  exit ( 0 );*/


  ps_WriteBBox ( 18, 12, 373, 286 );
  ps_OpenFile ( fn, 600 );

  DrawKnots ( N, u, 1 );

  ps_GSave ();
  ps_Write_Command ( "0 1100 translate" );
  ps_Set_Gray ( 0.68 );
  ps_Set_Line_Width ( 6.0 );
  DrawBSCurve ( n, N, u, d );
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  DrawCPolygon ( n, N, d );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i <= DD; i++ ) {
    t = u[n] + (float)i/(float)DD*(u[N-n]-u[n]);
    mbs_deBoorDerC2f ( n, N, u, (float*)d, t, &p, &d1 );
    ps_Fill_Circle ( p.x, p.y, 8.0 );
    psl_SetLine ( p.x, p.y, p.x+d1.x, p.y+d1.y, 0.0, 1.0 );
    psl_Draw ( 0.0, 0.97, 4.0 );
    psl_Arrow ( 1.0, true );
  }
  ps_GRestore ();

  ps_GSave ();
  ps_Write_Command ( "1500 1100 translate" );
  ps_Set_Gray ( 0.68 );
  ps_Set_Line_Width ( 6.0 );
  DrawBSCurve ( n, N, u, d );
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  DrawCPolygon ( n, N, d );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i <= DD; i++ ) {
    t = u[n] + (float)i/(float)DD*(u[N-n]-u[n]);
    mbs_deBoorDer2C2f ( n, N, u, d, t, &p, &d1, &d2 );
    ps_Fill_Circle ( p.x, p.y, 8.0 );
    psl_SetLine ( p.x, p.y, p.x+d2.x, p.y+d2.y, 0.0, 1.0 );
    psl_Draw ( 0.0, 0.97, 4.0 );
    psl_Arrow ( 1.0, true );
  }
  ps_GRestore ();

  ps_Set_Gray ( 0.68 );
  ps_Set_Line_Width ( 6.0 );
  DrawBSCurve ( n, N, u, d );
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  DrawCPolygon ( n, N, d );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i <= DD; i++ ) {
    t = u[n] + (float)i/(float)DD*(u[N-n]-u[n]);
    mbs_deBoorDer3C2f ( n, N, u, d, t, &p, &d1, &d2, &d3 );
    ps_Fill_Circle ( p.x, p.y, 8.0 );
    psl_SetLine ( p.x, p.y, p.x+d3.x, p.y+d3.y, 0.0, 1.0 );
    psl_Draw ( 0.0, 0.97, 4.0 );
    psl_Arrow ( 1.0, true );
  }

  ps_CloseFile ();
  printf ( "%s\n", fn );
  pkv_DestroyScratchMem ();
  exit ( 0 );
#undef DD
} /*main*/

