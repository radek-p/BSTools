
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

char fn[] = "intbsc.ps";


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

/* ///////////////////////////////////////////////////////////////////////// */

#define DEGREE    3
#define LASTIKNOT 10
float v[LASTIKNOT+1] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

point2f  p[LASTIKNOT+1];
vector2f da, db, dda, ddb;

float   u[100];
point2f d[100];
int lastknot;

static void SetupData ( void )
{
  int i;
  float a;

  for ( i = 0; i <= LASTIKNOT; i++ ) {
    a = (float)((((float)i+0.5)/(float)(LASTIKNOT+1))*2.0*PI);
    SetPoint2f ( &p[i], (float)(400.0*cos(a)+450.0), (float)(400.0*sin(a)+450.0) );
  }
  SetVector2f ( &da, 300.0 , 50.0 );
  SetVector2f ( &db, -300.0, 50.0 );
  SetVector2f ( &dda, 0.0, 0.0 );
  SetVector2f ( &ddb, 0.0, 0.0 );
} /*SetupData*/

static void DrawIPoints ( boolean ders )
{
  int i;

  for ( i = 0; i <= LASTIKNOT; i++ )
    ps_Fill_Circle ( p[i].x, p[i].y, 8.0 );

  if ( ders ) {
    psl_SetLine ( p[0].x, p[0].y, p[0].x+da.x, p[0].y+da.y, 0.0, 1.0 );
    psl_Draw ( 0.0, 0.9, 6.0 );
    psl_Arrow ( 1.0, true );
    psl_SetLine ( p[LASTIKNOT].x, p[LASTIKNOT].y,
                  p[LASTIKNOT].x+db.x, p[LASTIKNOT].y+db.y, 0.0, 1.0 );
    psl_Draw ( 0.0, 0.9, 6.0 );
    psl_Arrow ( 1.0, true );
  }
} /*DrawIPoints*/

static void Exp ( void )
{
  ps_WriteBBox ( 1, 3, 246, 346 );
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();

  SetupData ();
  ps_GSave ();
  ps_Write_Command ( "0 2000 translate" );
  mbs_multiBSCubicInterpf ( LASTIKNOT, v, 1, 2, 0, (float*)p, 0,
                            BS3_BC_FIRST_DER, (float*)&da,
                            BS3_BC_FIRST_DER, (float*)&db,
                            &lastknot, u, 0, (float*)d );
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
  mbs_multiBSCubicInterpf ( LASTIKNOT, v, 1, 2, 0, (float*)p, 0,
                            BS3_BC_BESSEL, NULL, BS3_BC_BESSEL, NULL,
                            &lastknot, u, 0, (float*)d );
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
  mbs_multiBSCubicInterpf ( LASTIKNOT, v, 1, 2, 0, (float*)p, 0,
                            BS3_BC_SECOND_DER, (float*)&da,
                            BS3_BC_SECOND_DER, (float*)&db,
                            &lastknot, u, 0, (float*)d );
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
  mbs_multiBSCubicInterpf ( LASTIKNOT, v, 1, 2, 0, (float*)p, 0,
                            BS3_BC_SECOND_DER0, NULL,
                            BS3_BC_SECOND_DER0, NULL,
                            &lastknot, u, 0, (float*)d );
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
  mbs_multiBSCubicInterpf ( LASTIKNOT, v, 1, 2, 0, (float*)p, 0,
                            BS3_BC_NOT_A_KNOT, NULL,
                            BS3_BC_NOT_A_KNOT, NULL,
                            &lastknot, u, 0, (float*)d );
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
  mbs_multiBSCubicInterpf ( LASTIKNOT, v, 1, 2, 0, (float*)p, 0,
                            BS3_BC_THIRD_DER0, NULL,
                            BS3_BC_THIRD_DER0, NULL,
                            &lastknot, u, 0, (float*)d );
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

