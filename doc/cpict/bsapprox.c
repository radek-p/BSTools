
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h> 
#include <string.h>
#include <math.h>

#include "multibs.h"
#include "psout.h"

/* #define DEBUG_PRINT */

char fn[] = "bsapprox.ps";

#define LASTQPOINT 16
#define LASTKNOT 13
#define DEGREE 3

float v[LASTQPOINT+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
point2f qc[LASTQPOINT+1] =
  {{300,550},{150,700},{150,850},{300,1000},
   {450,1150},{600,1150},{750,1000},{900,1150},
   {1050,1150},{1200,1000},{1350,1000},{1500,1000},
   {1500,850},{1350,850},{1200,700},{1050,550},{900,400}};

float u[LASTKNOT+1] = {-1,0,0,0,2,5,8,10,12,14,16,16,16,17};
point2f bsc[LASTKNOT-DEGREE];

bandm_profile aprof[LASTKNOT-DEGREE+1],
              qprof[LASTKNOT-DEGREE+1], rprof[LASTKNOT-DEGREE+1];
float *amat, *qmat, *rmat;


static void DrawPoints ( int npoints, point2f *p )
{
  int i;

  ps_Set_Line_Width ( 1.4 );
  ps_Set_Gray ( 0.5 );
  ps_Draw_Polyline2f ( npoints, p );
  for ( i = 0; i < npoints; i++ )
    ps_Fill_Circle ( p[i].x, p[i].y, 10.0 );
} /*DrawPoints*/

static void DrawArc ( int degree, int lastknot, float *knots, point2f *cpoints,
                      float t0, float t1 )
{
#define DD 150
  int i;
  float t;
  point2f c[DD+1];

  for ( i = 0; i <= DD; i++ ) {
    t = t0 + (t1-t0)*(float)i/DD;
    mbs_deBoorC2f ( degree, lastknot, knots, cpoints, t, &c[i] );
  }
  ps_Draw_Polyline2f ( DD+1, c );
#undef DD
} /*DrawArc*/

static void DrawCPoly ( int degree, int lastknot, point2f *cpoints )
{
  int i;

  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  ps_Draw_Polyline2f ( lastknot-degree, cpoints );
  for ( i = 0; i < lastknot-degree; i++ )
    ps_Mark_Circle ( cpoints[i].x, cpoints[i].y );
} /*DrawCPoly*/

static void DrawKnots ( int degree, int lastknot, const float *knots,
                        int lastpknot, const float *pknots )
{
  int   i;
  float x, y;

  psl_SetLine ( 100.0, 300.0, 1600.0, 300.0,
                knots[degree], knots[lastknot-degree] );
  psl_Draw ( knots[degree], knots[lastknot-degree], 2.0 );
  for ( i = 1; i < lastknot; i++ ) {
    if ( knots[i] > knots[i-1] ) {
      psl_Tick ( knots[i] );
      psl_GetPointf ( knots[i], &x, &y );
    }
    y -= 40.0;
    ps_Fill_Circle ( x, y, 10.0 );
  }
  for ( i = 0; i <= lastpknot; i++ )
    psl_HTick ( pknots[i], true );
} /*DrawKnots*/

int main ()
{
  ps_WriteBBox ( 11, 21, 196, 140 );
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();
  pkv_InitScratchMem ( 65536 );

  ps_Set_Gray ( 0.0 );
  DrawKnots ( DEGREE, LASTKNOT, u, LASTQPOINT, v );
  DrawPoints ( 16, qc );
  mbs_ConstructApproxBSC2f ( DEGREE, LASTKNOT, u, LASTQPOINT, v, qc, bsc );

  ps_Set_Gray ( 0.68 );
  ps_Set_Line_Width ( 6.0 );
  DrawArc ( DEGREE, LASTKNOT, u, bsc, u[DEGREE], u[LASTKNOT-DEGREE] );
  DrawCPoly ( DEGREE, LASTKNOT, bsc );

  pkv_DestroyScratchMem ();
  ps_CloseFile ();
  printf ( "%s\n", fn );
  exit ( 0 );
} /*main*/

