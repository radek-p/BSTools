
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "pknum.h"
#include "multibs.h"
#include "psout.h"

char fn[] = "bscfeded.ps";

#define n   3
#define NN 10
double   u[NN+1] = {0.0, 1.0, 3.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
point2d d[NN-n] = {{100.0,100.0},{100.0,400.0},{300.0,700.0},{800.0,900.0},
                   {1300.0,700.0},{1500.0,400.0},{1500.0,100.0}};

static void DrawCPolygon ( int degree, int lastknot, const point2d *d )
{
  int i;

  ps_Draw_Polyline2d ( lastknot-degree, d );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < lastknot-degree; i++ )
    ps_Mark_Circle ( (float)d[i].x, (float)d[i].y );
} /*DrawCPolygon*/

static void DrawBSCurve ( int degree, int lastknot, const double *u,
                          const point2d *d )
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

int main ()
{
  double   ui[50];
  point2d di[50];
  int     ni, NNi;

  pkv_InitScratchMem ( 65536 );
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();

  ps_Set_Gray ( 0.6 );
  ps_Set_Line_Width ( 8.0 );
  DrawCPolygon ( n, NN, d );
  DrawBSCurve ( n, NN, u, d );

  mbs_BSDegElevC2d ( n, NN, u, d, 1, &ni, &NNi, ui, di, false );
  printf ( "ni = %d, NNi = %d\n", ni, NNi );
  WriteArrayd ( "ui", NNi+1, ui );
  WriteArrayd ( "di", 2*(NNi-ni), &di[0].x );

  ps_Set_Line_Width ( 2.0 );
  ps_Write_Command ( "1 0 0 setrgbcolor" );
  DrawCPolygon ( ni, NNi, di );
  ps_Write_Command ( "0 1 0 setrgbcolor" );
  DrawBSCurve ( ni, NNi, ui, di );

  ps_CloseFile ();
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

