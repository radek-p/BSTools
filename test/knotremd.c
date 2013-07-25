
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"
#include "psout.h"

char fn[] = "knotremd.ps";

#define n 3
#define NN 12
int NNN = NN;

double kn[NN+1] = {-0.5,0.0,0.0,0.0,0.5,0.5,0.5,0.5,0.5,1.0,1.0,1.0,1.5};
point2d cp[NN-n] = {{1100.0,100.0},{1100.0,250.0},{1300.0,500.0},
                    {1800.0,700.0},{1900.0,400.0},{2200.0,400.0},
                    {2200.0,100.0},{2300.0,300.0},{2500.0,500.0}};

static void DrawArc ( double t0, double t1 )
{
#define DD 50
  int i;
  double t;
  point2d c[DD+1];

  for ( i = 0; i <= DD; i++ ) {
    t = t0 + (t1-t0)*(double)i/DD;
    mbs_deBoorC2d ( n, NN, kn, cp, t, &c[i] );
  }
  ps_Draw_Polyline2d ( DD+1, c );
#undef DD
} /*DrawArc*/

static void DrawKnots ( char c )
{
  int i;
  double w, x, y;

  if ( c ) {
    ps_Set_Gray ( 0.0 );  w = 2.0;
  }
  else {
    ps_Set_Gray ( 0.6 );  w = 6.0;
  }
  psl_SetLine ( 50.0, 300.0, 950.0, 300.0, (float)kn[1], (float)kn[NNN-1] );
  psl_Draw ( (float)kn[1], (float)kn[NNN-1], (float)w );
  for ( i = 1; i < NNN; i++ ) {
    if ( kn[i] > kn[i-1] ) {
      psl_Tick ( (float)kn[i] );
      psl_GetPointd ( kn[i], &x, &y );
    }
    y -= 40;
    ps_Fill_Circle ( (float)x, (float)y, 10 );
  }
} /*DrawKnots*/

static void ShowIt ( int dy )
{
  char s[30];
  int i;

  ps_GSave ();
  sprintf ( s, "0 %d translate", dy );
  ps_Write_Command ( s );
  DrawKnots ( 0 );
  ps_Set_Gray ( 0.6 );
  ps_Set_Line_Width ( 10.0 );
  DrawArc ( 0.0, 0.49999 );
  DrawArc ( 0.5, 1.0 );
  ps_Set_Line_Width ( 8.0 );
  ps_Draw_Polyline2d ( NNN-n, cp );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < NNN-n; i++ )
    ps_Mark_Circle ( (float)cp[i].x, (float)cp[i].y );

  mbs_KnotRemoveC2d ( n, &NNN, kn, cp, 5 );
  DrawKnots ( 1 );
  ps_Set_Line_Width ( 4.0 );
  DrawArc ( 0.0, 0.49999 );
  DrawArc ( 0.5, 1.0 );
  ps_Set_Line_Width ( 2.0 );
  ps_Draw_Polyline2d ( NNN-n, cp );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < NNN-n; i++ )
    ps_Mark_Circle ( (float)cp[i].x, (float)cp[i].y );
  ps_GRestore ();
} /*ShowIt*/

int main ()
{
  pkv_InitScratchMem ( 65536 );
  ps_WriteBBox ( 4, 12, 301, 302 );
  ps_OpenFile ( fn, 600 );

  ShowIt ( 1800 );
  ShowIt ( 1200 );
  ShowIt ( 600 );
  ShowIt ( 0 );

  ps_CloseFile ();
  pkv_DestroyScratchMem ( );
  exit ( 0 );
} /*main*/

