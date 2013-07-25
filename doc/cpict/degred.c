
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"
#include "psout.h"

char fn1[] = "degred1.ps";
char fn2[] = "degred2.ps";

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
                         float u0, float uN,
                         float ax, float bx, float ay, float by, char dom )
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
  psl_Arrow ( 1.35, true );
  psl_SetLine ( bx, by, bx+ax, by, u0, uN );
  psl_Draw ( u0, (float)(uN+0.3), 2.0 );
  psl_Draw ( knots[degree], knots[lastknot-degree], 6.0 );
  psl_Arrow ( (float)(uN+0.35), true );
  t = (float)(u0-1.0e10);
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
    if ( dom ) {
      ps_Set_Gray ( 0.0 );
      ps_Write_Command ( "0 setlinecap" );
      psl_SetLine ( bx, (float)(by-190.0-25*(i+1)), bx+ax,
                    (float)(by-190.0-25*(i+1)), u0, uN );
      psl_Draw ( knots[i], knots[i+degree+1], 6.0 );
      ps_Write_Command ( "1 setlinecap" );
      for ( j = i; j <= i+degree+1; j++ )
        psl_Tick ( knots[j] );
    }
    ps_Set_Gray ( 0.68 );
    ps_Set_Line_Width ( 6.0 );
    for ( j = 0; j <= DD; j++ ) {
      t = knots[i] + (float)j/(float)DD*(knots[i+degree+1]-knots[i]);
      mbs_deBoorBasisf ( degree, lastknot, knots, t, &fnz, &nnz, f );
      c[j].x = bx+ax*(t-u0)/(uN-u0);
      c[j].y = by;
      if ( i >= fnz && i < fnz+nnz )
        c[j].y += ay*f[i-fnz];
    }
    ps_Draw_Polyline2f ( DD+1, c );
  }
  ps_Set_Gray ( 0.0 );
  psl_SetLine ( bx, by, bx+ax, by, u0, uN );
  psl_Draw ( knots[degree], knots[lastknot-degree], 6.0 );
#undef DD
} /*DrawGraphs*/

/* this one is closed */
#define deg  5
#define LKN 16
float knots[17] = {-0.668137,-0.343109,-0.343109,-0.191582,0.111930,
  0.111930,0.331863,0.656891,0.656891,0.808418,1.111930,1.111930,
  1.331863,1.656891,1.656891,1.808418,2.111930};
point2f cpoints[11] =
  {{0.639733,0.771523},{0.861817,0.278146},{0.570125,-0.642384},
   {-0.132587,-0.615894},{-0.997718,-0.205298},{-0.507146,0.857616},
   {0.639733,0.771523},{0.861817,0.278146},{0.570125,-0.642384},
   {-0.132587,-0.615894},{-0.997718,-0.205298}};


float outknots[50];
point2f outcpoints[50];
int outN, outdeg;

int main ( void )
{
  int i;

  pkv_InitScratchMem ( 65536 );

  for ( i = 0; i < LKN-deg; i++ ) {
    cpoints[i].x = (float)(500.0*(cpoints[i].x+1.25));
    cpoints[i].y = (float)(500.0*(cpoints[i].y+1.2));
  }

  ps_WriteBBox ( 13, 16, 331, 147 );
  ps_OpenFile ( fn1, 600 );

  mbs_multiBSDegRedf ( 1, 2, deg, LKN, knots, 0, (float*)cpoints, 1,
                       &outdeg, &outN, outknots, 0, (float*)outcpoints );

  ps_Set_Gray ( 0.72 );
  ps_Set_Line_Width ( 8.0 );
  DrawBSCurve ( deg, LKN, knots, cpoints );
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  DrawBSCurve ( outdeg, outN, outknots, outcpoints );
  ps_Set_Line_Width ( 4.0 );
  ps_Set_Gray ( 0.6 );
  DrawBSPolyline ( deg, LKN, cpoints, 0 );
  ps_Set_Line_Width ( 2.0 );
  ps_Set_Gray ( 0.0 );
  DrawBSPolyline ( outdeg, outN, outcpoints, 1 );
  DrawGraphs ( deg, LKN, knots, knots[0], knots[LKN],
               1200, 1400, 300, 800, false );
  DrawGraphs ( outdeg, outN, outknots, knots[0], knots[LKN],
               1200, 1400, 300, 200, false );
  ps_CloseFile ();
  printf ( "%s\n", fn1 );

  ps_WriteBBox ( 13, 16, 331, 145 );
  ps_OpenFile ( fn2, 600 );
  mbs_multiBSDegRedClosedf ( 1, 2, deg, LKN, knots, 0, (float*)cpoints, 1,
                       &outdeg, &outN, outknots, 0, (float*)outcpoints );

  ps_Set_Gray ( 0.72 );
  ps_Set_Line_Width ( 8.0 );
  DrawBSCurve ( deg, LKN, knots, cpoints );
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  DrawBSCurve ( outdeg, outN, outknots, outcpoints );
  ps_Set_Line_Width ( 4.0 );
  ps_Set_Gray ( 0.6 );
  DrawBSPolyline ( deg, LKN, cpoints, 0 );
  ps_Set_Line_Width ( 2.0 );
  ps_Set_Gray ( 0.0 );
  DrawBSPolyline ( outdeg, outN, outcpoints, 1 );
  DrawGraphs ( deg, LKN, knots, knots[0], knots[LKN],
               1200, 1400, 300, 800, false );
  DrawGraphs ( outdeg, outN, outknots, knots[0], knots[LKN],
               1200, 1400, 300, 200, false );

  ps_CloseFile ();
  printf ( "%s\n", fn2 );
  printf ( "Scratch memory used: %d\n", (int)pkv_MaxScratchTaken () );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

