
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

char fn[] = "bsder123d.ps";

static void DrawCPolygon ( int degree, int lastknot, const point2d *d )
{
  int i;

  ps_Draw_Polyline2d ( lastknot-degree, d );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < lastknot-degree; i++ )
    ps_Mark_Circle ( (float)d[i].x, (float)d[i].y );
} /*DrawCPolygon*/

static void DrawBSCurve ( int degree, int lastknot, const double *u, const point2d *d )
{
#define DD 200
  int     i;
  double   t;
  point2d *c;

  c = pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  for ( i = 0; i <= DD; i++ ) {
    t = u[degree] + (double)i/(double)DD*(u[lastknot-degree]-u[degree]);
    mbs_deBoorC2d ( degree, lastknot, u, d, t, &c[i] );
  }
  ps_Draw_Polyline2d ( DD+1, c );
  pkv_FreeScratchMem ( (DD+1)*sizeof(point2d) );
#undef DD
} /*DrawBSCurve*/

static void DrawKnots ( int NNN, double *kn, char c )
{
  int i;
  double w, x, y;

  if ( c ) {
    ps_Set_Gray ( 0.0 );  w = 2.0;
  }
  else {
    ps_Set_Gray ( 0.6 );  w = 6.0;
  }
  psl_SetLine ( 1600.0, 350.0, 3000.0, 350.0, (float)kn[1], (float)kn[NNN-1] );
  psl_Draw ( (float)kn[1], (float)kn[NNN-1], (float)w );
  for ( i = 1; i < NNN; i++ ) {
    if ( i == 1 || kn[i] > kn[i-1] ) {
      psl_Tick ( (float)kn[i] );
      psl_GetPointd ( kn[i], &x, &y );
    }
    y -= 40;
    ps_Fill_Circle ( (float)x, (float)y, 10 );
  }
} /*DrawKnots*/

/*
static void Test1 ( double gr, int n, double *d )
{
  int i;

  ps_Set_Gray ( gr );
  ps_Set_Line_Width ( 3.0 );
  if ( n )
    ps_Draw_Polyline2d ( n+1, (point2d*)d );
  for ( i = 0; i <= n; i++ )
    ps_Fill_Circle ( d[i+i], d[i+i+1], 8.0 );
}*/ /*Test1*/

/* ///////////////////////////////////////////////////////////////////////// */
#define N 15
#define n  5
double u[N+1] = {0.0,0.0,0.0,0.0,0.0,0.0,
                3.5,3.5,3.5,3.5,
                7.0,7.0,7.0,7.0,7.0,7.0};

point2d d[N-n] =
  {{640.0,400.0},{320.0,400.0},{160.0,720.0},{640.0,1200.0},{1120.0,1200.0},
   {1440.0,880.0},{1120.0,720.0},{1120.0,400.0},{1440.0,400.0},{1600.0,560.0}};

point2d  p;
vector2d d1, d2, d3;

int main ()
{
#define DD 60
  int   i;
  double t;

  pkv_InitScratchMem ( 65536 );

  mbs_multideBoorDer3d ( n, N, u, 1, 2, 0, (double*)d,
                       3.5, (double*)&p, (double*)&d1, (double*)&d2, (double*)&d3 );
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
    t = u[n] + (double)i/(double)DD*(u[N-n]-u[n]);
    mbs_deBoorDerC2d ( n, N, u, d, t, &p, &d1 );
    ps_Fill_Circle ( (float)p.x, (float)p.y, 8.0 );
    psl_SetLine ( (float)p.x, (float)p.y, (float)(p.x+d1.x),
                  (float)(p.y+d1.y), 0.0, 1.0 );
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
  ps_Set_RGB ( 1.0, 0.0, 0.0 );
  for ( i = 0; i <= DD; i++ ) {
    t = u[n] + (double)i/(double)DD*(u[N-n]-u[n]);
    mbs_deBoorDer2C2d ( n, N, u, (double*)d, t, &p, &d1, &d2 );
    ps_Fill_Circle ( (float)p.x, (float)p.y, 8.0 );
    psl_SetLine ( (float)p.x, (float)p.y, (float)(p.x+d1.x),
                  (float)(p.y+d1.y), 0.0, 1.0 );
    psl_Draw ( 0.0, 0.97, 4.0 );
    psl_Arrow ( 1.0, true );
  }
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i <= DD; i++ ) {
    t = u[n] + (double)i/(double)DD*(u[N-n]-u[n]);
    mbs_deBoorDer2C2d ( n, N, u, d, t, &p, &d1, &d2 );
    ps_Fill_Circle ( (float)p.x, (float)p.y, 8.0 );
    psl_SetLine ( (float)p.x, (float)p.y, (float)(p.x+d2.x),
                  (float)(p.y+d2.y), 0.0, 1.0 );
    psl_Draw ( 0.0, 0.97, 4.0 );
    psl_Arrow ( 1.0, true );
    psl_SetLine ( (float)p.x, (float)p.y, (float)(p.x+d1.x),
                  (float)(p.y+d1.y), 0.0, 1.0 );
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
  ps_Set_RGB ( 1.0, 0.0, 0.0 );
  for ( i = 0; i <= DD; i++ ) {
    t = u[n] + (double)i/(double)DD*(u[N-n]-u[n]);
    mbs_deBoorDer2C2d ( n, N, u, (double*)d, t, &p, &d1, &d2 );
    ps_Fill_Circle ( (float)p.x, (float)p.y, 8.0 );
    psl_SetLine ( (float)p.x, (float)p.y, (float)(p.x+d2.x),
                  (float)(p.y+d2.y), 0.0, 1.0 );
    psl_Draw ( 0.0, 0.97, 4.0 );
    psl_Arrow ( 1.0, true );
  }
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i <= DD; i++ ) {
    t = u[n] + (double)i/(double)DD*(u[N-n]-u[n]);
    mbs_deBoorDer3C2d ( n, N, u, d, t, &p, &d1, &d2, &d3 );
    ps_Fill_Circle ( (float)p.x, (float)p.y, 8.0 );
    psl_SetLine ( (float)p.x, (float)p.y, (float)(p.x+d3.x),
                  (float)(p.y+d3.y), 0.0, 1.0 );
    psl_Draw ( 0.0, 0.97, 4.0 );
    psl_Arrow ( 1.0, true );
    psl_SetLine ( (float)p.x, (float)p.y, (float)(p.x+d2.x),
                  (float)(p.y+d2.y), 0.0, 1.0 );
    psl_Draw ( 0.0, 0.97, 4.0 );
    psl_Arrow ( 1.0, true );
  }

  ps_CloseFile ();
  printf ( "%s\n", fn );
  pkv_DestroyScratchMem ();
  exit ( 0 );
#undef DD
} /*main*/

