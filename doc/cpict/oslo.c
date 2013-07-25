
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

char fn[] = "oslo.ps";

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

static void DrawGraphs ( int degree, int lastknot, const float *knots,
                         float ax, float bx, float ay, float by )
{
#define DD 50
  int     i, j, fnz, nnz;
  float   t, x, y, f[10];
  point2f c[DD+1];

  ps_Set_Gray ( 0.0 );
  psl_SetLine ( bx, by, bx, by+ay, 0.0, 1.0 );
  psl_Draw ( 0.0, 1.2, 2.0 );
  psl_Tick ( 0.0 );
  psl_Tick ( 1.0 );
  psl_Arrow ( 1.25, true );
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
  for ( i = 0; i < lastknot-degree; i++ ) {
    ps_Set_Gray ( 0.0 );
    ps_Write_Command ( "0 setlinecap" );
    psl_SetLine ( bx, (float)(by-150.0-25*(i+1)), bx+ax,
                 (float)(by-150.0-25*(i+1)), knots[0], knots[lastknot]);
    psl_Draw ( knots[i], knots[i+degree+1], 6.0 );
    ps_Write_Command ( "1 setlinecap" );
    for ( j = i; j <= i+degree+1; j++ )
      psl_Tick ( knots[j] );
    ps_Set_Gray ( 0.68 );
    ps_Set_Line_Width ( 6.0 );
    for ( j = 0; j <= DD; j++ ) {
      t = knots[i] + (float)j/(float)DD*(knots[i+degree+1]-knots[i]);
      mbs_deBoorBasisf ( degree, lastknot, knots, t, &fnz, &nnz, f );
      c[j].x = bx+ax*(t-knots[0])/(knots[lastknot]-knots[0]);
      c[j].y = by;
      if ( i >= fnz && i < fnz+nnz )
        c[j].y += ay*f[i-fnz];
    }
    ps_Draw_Polyline2f ( DD+1, c );
  }
#undef DD
} /*DrawGraphs*/

/* ///////////////////////////////////////////////////////////////////////// */

#define DEGREE    3
#define LASTUKNOT 14
#define LASTVKNOT 18
float u[LASTUKNOT+1] = {0.0, 1.0, 1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                        7.0, 8.0, 9.0, 9.0, 9.0, 10.0};
float v[LASTVKNOT+1] = {0.0, 1.0, 1.0, 1.0, 2.0, 3.0, 3.5, 4.0, 4.5,
                        5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 9.0, 9.0, 10.0};

point2f d[LASTUKNOT-DEGREE] =
  {{100.0,2100.0},{100.0,2400.0},{400.0,2500.0},{800.0,2400.0},{300.0,2200.0},
   {800.0,2100.0},{900.0,2600.0},{1200.0,2200.0},{1700.0,2800.0},
   {2550.0,2500.0},{2200.0,2100.0}};

point2f da[LASTVKNOT-DEGREE];

int main ()
{
  int           nc, as;
  float         *a;
  bandm_profile *prof;

  ps_WriteBBox ( 10, 8, 308, 338 );
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();

  pkv_InitScratchMem ( 65536 );
  nc = LASTUKNOT-DEGREE;
  prof = pkv_GetScratchMem ( (nc+1)*sizeof(bandm_profile) );

  if ( mbs_OsloKnotsCorrectf ( LASTUKNOT, u, LASTVKNOT, v ) ) {
    as = mbs_BuildOsloMatrixProfilef ( DEGREE, LASTUKNOT, u, LASTVKNOT, v,
                                       prof );

    a = pkv_GetScratchMem ( as*sizeof(float) );
    mbs_BuildOsloMatrixf ( DEGREE, LASTUKNOT, u, v, prof, a );
/*
    pkn_PrintBandmf ( nc, prof, a );
*/
    pkn_multiBandmMultVectorf ( LASTVKNOT-DEGREE, LASTUKNOT-DEGREE,
                                prof, a, 2, (float*)d, (float*)da );

    ps_Set_Gray ( 0.0 );
    ps_Set_Line_Width ( 8.0 );
    DrawBSCurve ( DEGREE, LASTUKNOT, u, d );
    ps_Set_Gray ( 0.72 );
    ps_Set_Line_Width ( 2.0 );
    DrawBSCurve ( DEGREE, LASTVKNOT, v, da );

    ps_Set_Line_Width ( 4.0 );
    ps_Set_Gray ( 0.0 );
    DrawBSPolyline ( DEGREE, LASTUKNOT, d, 0 );

    ps_Set_Line_Width ( 2.0 );
    ps_Set_Gray ( 0.68 );
    DrawBSPolyline ( DEGREE, LASTVKNOT, da, 1 );

    pkv_FreeScratchMem ( as*sizeof(float) );
  }
  DrawGraphs ( DEGREE, LASTUKNOT, u, 2200.0, 100.0, 300.0, 1600.0 );
  DrawGraphs ( DEGREE, LASTVKNOT, v, 2200.0, 100.0, 300.0, 600.0 );

  pkv_DestroyScratchMem ();
  ps_CloseFile ();
  printf ( "%s\n", fn );
  exit ( 0 );
} /*main*/

