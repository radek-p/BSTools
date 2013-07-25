
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

char fn[] = "skiplrd.ps";


/* ///////////////////////////////////////////////////////////////////////// */
static void DrawBSPolyline ( int degree, int lastknot, const point2d *d, int mk )
{
  int i;

  ps_Draw_Polyline2d ( lastknot-degree, d );

  switch ( mk ) {
case 0:
    for ( i = 0; i < lastknot-degree; i++ )
      ps_Mark_Circle ( (float)d[i].x, (float)d[i].y );
    break;
case 1:
    for ( i = 0; i < lastknot-degree; i++ )
      ps_Fill_Circle ( (float)d[i].x, (float)d[i].y, 8.0 );
    break;
default:
    ;
  }
} /*DrawBSPolyline*/

static void DrawBSCurve ( int degree, int lastknot,
                          const double *knots, const point2d *d )
{
#define DD 32
/*
  int i;
  point2d c[DD+1];
  double t;

  for ( i = 0; i <= DD; i++ ) {
    t = knots[degree]+(double)i/(double)DD*(knots[lastknot-degree]-knots[degree]);
    mbs_deBoorC2d ( degree, lastknot, knots, d, t, &c[i] );
  }
  ps_Draw_Polylined ( c, DD );
*/

  int     i, j, kpcs;
  point2d *ncp, c[DD+1];
  int     scratchsize;
  double   t;
  
  kpcs = mbs_NumKnotIntervalsd ( degree, lastknot, knots );
  ncp = (point2d*)pkv_GetScratchMem (
                    scratchsize = (degree+1)*kpcs*sizeof(point2d) );
  mbs_BSToBezC2d ( degree, lastknot, knots, d, &kpcs, NULL, NULL, ncp );
  for ( i = 0; i < kpcs; i++ ) {
    for ( j = 0; j <= DD; j++ ) {
      t = (double)j/(double)DD;
      mbs_BCHornerC2d ( degree, &ncp[i*(degree+1)], t, &c[j] );
    }
    ps_Draw_Polyline2d ( DD+1, c );
  }
  pkv_FreeScratchMem ( scratchsize );

#undef DD
} /*DrawBSCurve*/

static void DrawGraphs ( int degree, int lastknot, const double *knots,
                         double ax, double bx, double ay, double by, char dom )
{
#define DD 50
  int     i, j, fnz, nnz;
  double   t, x, y, f[10];
  point2d c[DD+1];

  ps_Set_Gray ( 0.0 );
  psl_SetLine ( (float)bx, (float)by, (float)bx, (float)(by+ay), 0.0, 1.0 );
  psl_Draw ( 0.0, 1.2, 2.0 );
  psl_Tick ( 0.0 );
  psl_Tick ( 1.0 );
  psl_Arrow ( 1.25, true );
  psl_SetLine ( (float)bx, (float)by, (float)(bx+ax), (float)by,
                (float)knots[0], (float)knots[lastknot] );
  psl_Draw ( (float)knots[0], (float)(knots[lastknot]+1.0), 2.0 );
  psl_Arrow ( (float)(knots[lastknot]+1.1), true );
  t = (float)(knots[0]-1.0e10);
  for ( i = 0; i <= lastknot; i++ ) {
    psl_Tick ( (float)knots[i] );
    if ( knots[i] > t ) {
      psl_GetPointd ( knots[i], &x, &y );
    }
    y -= 40.0;
    ps_Fill_Circle ( (float)x, (float)y, 10.0 );
    t = knots[i];
  }
  for ( i = 0; i < lastknot-degree; i++ ) {
    if ( dom ) {
      ps_Set_Gray ( 0.0 );
      ps_Write_Command ( "0 setlinecap" );
      psl_SetLine ( (float)bx, (float)(by-190.0-25*(i+1)), (float)(bx+ax),
                    (float)(by-190.0-25*(i+1)), (float)knots[0], (float)knots[lastknot]);
      psl_Draw ( (float)knots[i], (float)knots[i+degree+1], 6.0 );
      ps_Write_Command ( "1 setlinecap" );
      for ( j = i; j <= i+degree+1; j++ )
        psl_Tick ( (float)knots[j] );
    }
    ps_Set_Gray ( 0.68 );
    ps_Set_Line_Width ( 6.0 );
    for ( j = 0; j <= DD; j++ ) {
      t = knots[i] + (double)j/(double)DD*(knots[i+degree+1]-knots[i]);
      mbs_deBoorBasisd ( degree, lastknot, knots, t, &fnz, &nnz, f );
      c[j].x = bx+ax*(t-knots[0])/(knots[lastknot]-knots[0]);
      c[j].y = by;
      if ( i >= fnz && i < fnz+nnz )
        c[j].y += ay*f[i-fnz];
    }
    ps_Draw_Polyline2d ( DD+1, c );
  }
#undef DD
} /*DrawGraphs*/

/* ///////////////////////////////////////////////////////////////////////// */

#define DEGREE    3
#define LASTUKNOT 10
double u[LASTUKNOT+1] = {-1.0, 0.0, 0.0, 1.0, 1.0, 1.5, 2.0, 2.0, 3.0, 3.0, 4.0};

point2d d[LASTUKNOT-DEGREE] =
  {{100.0,2100.0},{100.0,2400.0},{400.0,2500.0},{500.0,2600.0},
   {800.0,2400.0},{300.0,2200.0},{800.0,2100.0}};

static void Exp ()
{
  ps_WriteBBox ( 10, 100, 308, 338 );
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();

  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 8.0 );
  DrawBSCurve ( DEGREE, LASTUKNOT, u, d );

  ps_Set_Line_Width ( 4.0 );
  ps_Set_Gray ( 0.0 );
  DrawBSPolyline ( DEGREE, LASTUKNOT, d, 0 );

  DrawGraphs ( DEGREE, LASTUKNOT, u, 2200, 100, 300, 1600, false );

  ps_CloseFile ();
  printf ( "%s\n", fn );
} /*Exp*/

int main ()
{
  pkv_InitScratchMem ( 65536 );
  Exp ();
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

