
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "multibs.h"
#include "psout.h"

char fn[] = "addspld.ps";

#define n1   3
#define NN1  9
#define n2   4
#define NN2 13

double kn1[NN1+1] =
  { 0.0, 0.0, 0.0, 0.0, 1.0/3.0, 2.0/3.0, 1.0, 1.0, 1.0, 1.0 };

double d1[NN1-n1] = {1.0,1.0,2.0,1.0,2.0,2.0};

double kn2[NN2+1] =
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 2.0/3.0,
    1.0, 1.0, 1.0, 1.0, 1.0 };

double d2[NN2-n2] = {0.5,0.0,0.5,1.0,1.5,1.0,2.0,2.5,2.5};

int deg, slkn;
double sumknots[50];
double scp[50];

double xmi, xma, ax, bx, ay, by;

static void SetupGraph ( double xmin, double xmax, double ymin, double ymax,
                  double ximin, double ximax, double etamin, double etamax )
{
  xmi = xmin;
  xma = xmax;
  ax = (ximin-ximax)/(xmin-xmax);
  bx = ximin - ax*xmin;
  ay = (etamin-etamax)/(ymin-ymax);
  by = etamin - ay*ymin;
  psl_SetLine ( (float)ximin, (float)etamin, (float)ximax, (float)etamin, 0.0, 1.0 );
  psl_Draw ( 0.0, 1.1, 2.0 );
  psl_Arrow ( 1.15, true );
  psl_SetLine ( (float)ximin, (float)etamin, (float)ximin, (float)etamax, 0.0, 1.0 );
  psl_Draw ( 0.0, 1.1, 2.0 );
  psl_Arrow ( 1.15, true );
} /*SetupGraph*/

static void DrawGraph ( int n, int NN, const double *knots, const double *coeff )
{
#define DD 50
  void    *sp;
  int     i;
  double   u, f;
  point2d *cc;

  sp = pkv_GetScratchMemTop ();
  cc = pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  if ( !cc )
    exit ( 1 );
  for ( i = 0; i <= DD; i++ ) {
    u = xmi + (double)i/(double)DD*(xma-xmi);
    mbs_deBoorC1d ( n, NN, knots, coeff, u, &f );
    SetPoint2d ( &cc[i], ax*u+bx, ay*f+by );
  }
  ps_Draw_Polyline2d ( DD+1, cc );
  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawGraph*/

static void DrawPolyline ( int n, int NN, const double *knots, const double *coeff )
{
  void    *sp;
  int     i, j;
  point2d *cc;
  double   xi;

  sp = pkv_GetScratchMemTop ();
  cc = pkv_GetScratchMem ( (NN-n)*sizeof(point2d) );
  if ( !cc )
    exit ( 1 );
  for ( i = 0; i < NN-n; i++ ) {
    xi = knots[i+1];
    for ( j = 2; j <= n; j++ )
      xi += knots[i+j];
    xi /= (double)n;
    SetPoint2d ( &cc[i], ax*xi+bx, ay*coeff[i]+by );
  }
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  ps_Draw_Polyline2d ( NN-n, cc );
  for ( i = 0; i < NN-n; i++ )
    ps_Mark_Circle ( (float)cc[i].x, (float)cc[i].y );
  pkv_SetScratchMemTop ( sp );
} /*DrawPolyline*/

static void DrawChkGraph ( void )
{
#define DD 50
  void    *sp;
  int     i;
  double   u, f1, f2;
  point2d *cc;

  sp = pkv_GetScratchMemTop ();
  cc = pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  if ( !cc )
    exit ( 1 );
  for ( i = 0; i <= DD; i++ ) {
    u = xmi + (double)i/(double)DD*(xma-xmi);
    mbs_deBoorC1d ( n1, NN1, kn1, d1, u, &f1 );
    mbs_deBoorC1d ( n2, NN2, kn2, d2, u, &f2 );
    SetPoint2d ( &cc[i], ax*u+bx, ay*(f1+f2)+by );
  }
  ps_Draw_Polyline2d ( DD+1, cc );
  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawChkGraph*/

int main ()
{
  pkv_InitScratchMem ( 65536 );
  ps_OpenFile ( fn, 600 );
  ps_Write_Command ( "1 setlinecap" );
  memset ( sumknots, 0, 50*sizeof(double) );

  mbs_multiAddBSCurvesd ( 1, 1, n1, NN1, kn1, 0, d1, n2, NN2, kn2, 0, d2,
                          &deg, &slkn, sumknots, 0, scp );

  SetupGraph ( 0.0, 1.0, 0.0, 4.5, 100.0, 1300.0, 100.0, 900.0 );

  DrawPolyline ( n1, NN1, kn1, d1 );
  ps_Set_Gray ( 0.68 );
  ps_Set_Line_Width ( 4.0 );
  DrawGraph ( n1, NN1, kn1, d1 );

  DrawPolyline ( n2, NN2, kn2, d2 );
  ps_Set_Gray ( 0.44 );
  ps_Set_Line_Width ( 4.0 );
  DrawGraph ( n2, NN2, kn2, d2 );

  DrawPolyline ( deg, slkn, sumknots, scp );
  ps_Set_Gray ( 0.72 );
  ps_Set_Line_Width ( 8.0 );
  DrawGraph ( deg, slkn, sumknots, scp );
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  DrawChkGraph ();

  ps_CloseFile ();
  printf ( "%s\n", fn );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main */

