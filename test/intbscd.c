
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

char fn[] = "intbscd.ps";


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
#define DD 200
  int i;
  point2d c[DD+1];
  double t;

  for ( i = 0; i <= DD; i++ ) {
    t = knots[degree]+(double)i/(double)DD*(knots[lastknot-degree]-knots[degree]);
    mbs_deBoorC2d ( degree, lastknot, knots, d, t, &c[i] );
  }
  ps_Draw_Polyline2d ( DD+1, c );
#undef DD
} /*DrawBSCurve*/

/* ///////////////////////////////////////////////////////////////////////// */

#define DEGREE    3
#define LASTIKNOT 10
double v[LASTIKNOT+1] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

point2d  p[LASTIKNOT+1];
vector2d da, db, dda, ddb;

double   u[100];
point2d d[100];
int lastknot;

static void SetupData ( void )
{
  int i;
  double a;

  for ( i = 0; i <= LASTIKNOT; i++ ) {
    a = (((double)i+0.5)/(double)(LASTIKNOT+1))*2.0*PI;
    SetPoint2d ( &p[i], 400.0*cos(a)+450, 400.0*sin(a)+450 );
  }
  SetVector2d ( &da, 300.0 , 50.0 );
  SetVector2d ( &db, -300.0, 50.0 );
  SetVector2d ( &dda, 0.0, 0.0 );
  SetVector2d ( &ddb, 0.0, 0.0 );
} /*SetupData*/

static void DrawIPoints ( boolean ders )
{
  int i;

  for ( i = 0; i <= LASTIKNOT; i++ )
    ps_Fill_Circle ( (float)p[i].x, (float)p[i].y, 8.0 );

  if ( ders ) {
    psl_SetLine ( (float)p[0].x, (float)p[0].y, (float)(p[0].x+da.x),
                  (float)(p[0].y+da.y), 0.0, 1.0 );
    psl_Draw ( 0.0, 0.9, 6.0 );
    psl_Arrow ( 1.0, true );
    psl_SetLine ( (float)p[LASTIKNOT].x, (float)p[LASTIKNOT].y,
                  (float)(p[LASTIKNOT].x+db.x), (float)(p[LASTIKNOT].y+db.y), 0.0, 1.0 );
    psl_Draw ( 0.0, 0.9, 6.0 );
    psl_Arrow ( 1.0, true );
  }
} /*DrawIPoints*/

static void Exp ()
{
  ps_WriteBBox ( 1, 3, 280, 346 );
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();

  SetupData ();
  ps_GSave ();
  ps_Write_Command ( "0 2000 translate" );
  mbs_multiBSCubicInterpd ( LASTIKNOT, v, 1, 2, 0, (double*)p, 0,
                            BS3_BC_FIRST_DER, (double*)&da,
                            BS3_BC_FIRST_DER, (double*)&db,
                            &lastknot, u, 0, (double*)d );
  ps_Set_Gray ( 0.68 );
  ps_Set_Line_Width ( 6.0 );
  DrawBSCurve ( 3, lastknot, u, d );
  ps_Set_Gray ( 0.0 );
  DrawIPoints ( true );
  ps_Set_Line_Width ( 2.0 );
  DrawBSPolyline ( 3, lastknot, d, 0 );
  ps_GRestore ();

  ps_GSave ();
  ps_Write_Command ( "1200 2000 translate" );
  mbs_multiBSCubicInterpd ( LASTIKNOT, v, 1, 2, 0, (double*)p, 0,
                            BS3_BC_BESSEL, NULL, BS3_BC_BESSEL, NULL,
                            &lastknot, u, 0, (double*)d );
  ps_Set_Gray ( 0.68 );
  ps_Set_Line_Width ( 6.0 );
  DrawBSCurve ( 3, lastknot, u, d );
  ps_Set_Gray ( 0.0 );
  DrawIPoints ( false );
  ps_Set_Line_Width ( 2.0 );
  DrawBSPolyline ( 3, lastknot, d, 0 );
  ps_GRestore ();
  
  ps_GSave ();
  ps_Write_Command ( "0 1000 translate" );
  mbs_multiBSCubicInterpd ( LASTIKNOT, v, 1, 2, 0, (double*)p, 0,
                            BS3_BC_SECOND_DER, (double*)&da,
                            BS3_BC_SECOND_DER, (double*)&db,
                            &lastknot, u, 0, (double*)d );
  ps_Set_Gray ( 0.68 );
  ps_Set_Line_Width ( 6.0 );
  DrawBSCurve ( 3, lastknot, u, d );
  ps_Set_Gray ( 0.0 );
  DrawIPoints ( true );
  ps_Set_Line_Width ( 2.0 );
  DrawBSPolyline ( 3, lastknot, d, 0 );
  ps_GRestore ();

  ps_GSave ();
  ps_Write_Command ( "1200 1000 translate" );
  mbs_multiBSCubicInterpd ( LASTIKNOT, v, 1, 2, 0, (double*)p, 0,
                            BS3_BC_SECOND_DER0, NULL,
                            BS3_BC_SECOND_DER0, NULL,
                            &lastknot, u, 0, (double*)d );
  ps_Set_Gray ( 0.68 );
  ps_Set_Line_Width ( 6.0 );
  DrawBSCurve ( 3, lastknot, u, d );
  ps_Set_Gray ( 0.0 );
  DrawIPoints ( false );
  ps_Set_Line_Width ( 2.0 );
  DrawBSPolyline ( 3, lastknot, d, 0 );
  ps_GRestore ();

  ps_GSave ();
  ps_Write_Command ( "0 0 translate" );
  mbs_multiBSCubicInterpd ( LASTIKNOT, v, 1, 2, 0, (double*)p, 0,
                            BS3_BC_NOT_A_KNOT, NULL,
                            BS3_BC_NOT_A_KNOT, NULL,
                            &lastknot, u, 0, (double*)d );
  ps_Set_Gray ( 0.68 );
  ps_Set_Line_Width ( 6.0 );
  DrawBSCurve ( 3, lastknot, u, d );
  ps_Set_Gray ( 0.0 );
  DrawIPoints ( false );
  ps_Set_Line_Width ( 2.0 );
  DrawBSPolyline ( 3, lastknot, d, 0 );
  ps_GRestore ();

  ps_GSave ();
  ps_Write_Command ( "1200 0 translate" );
  mbs_multiBSCubicInterpd ( LASTIKNOT, v, 1, 2, 0, (double*)p, 0,
                            BS3_BC_THIRD_DER0, NULL,
                            BS3_BC_THIRD_DER0, NULL,
                            &lastknot, u, 0, (double*)d );
  ps_Set_Gray ( 0.68 );
  ps_Set_Line_Width ( 6.0 );
  DrawBSCurve ( 3, lastknot, u, d );
  ps_Set_Gray ( 0.0 );
  DrawIPoints ( false );
  ps_Set_Line_Width ( 2.0 );
  DrawBSPolyline ( 3, lastknot, d, 0 );
  ps_GRestore ();

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

