
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

char fn[] = "bsderd.ps";

#define n 3
#define NN 12
int NNN = NN;

double kn[NN+1] = {-0.5,0.0,0.0,0.0,2.4,2.4,2.4,2.4,2.4,6.0,6.0,6.0,9.0};
point2d cp[NN-n] = {{1100.0,100.0},{1100.0,250.0},{1300.0,500.0},
                    {1800.0,700.0},{1900.0,400.0},{2200.0,400.0},
                    {2200.0,100.0},{2300.0,300.0},{2500.0,500.0}};
vector2d dcp[NN-n-1];
double dkn[NN];

static void DrawArc ( int nn, int NNN, double *kn, point2d *cp,
                      double t0, double t1 )
{
#define DD 50
  int i;
  double t;
  point2d c[DD+1];

  for ( i = 0; i <= DD; i++ ) {
    t = t0 + (t1-t0)*(double)i/DD;
    mbs_deBoorC2d ( nn, NNN, kn, cp, t, &c[i] );
  }
  ps_Draw_Polyline2d ( DD+1, c );
#undef DD
} /*DrawArc*/

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
  psl_SetLine ( 50.0, 300.0, 950.0, 300.0, (float)kn[1], (float)kn[NNN-1] );
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

static void ShowIt ( int dy )
{
  char s[30];
  int i;

  ps_GSave ();
  sprintf ( s, "0 %d translate", dy );
  ps_Write_Command ( s );
  DrawKnots ( NN, kn, 0 );
  ps_Set_Line_Width ( 6.0 );
  ps_Set_Gray ( 0.6 );
  DrawArc ( n, NN, kn, cp, 0.0, 2.3999 );
  DrawArc ( n, NN, kn, cp, 2.4, 6.0 );
  ps_Set_Line_Width ( 3.0 );
  ps_Set_Gray ( 0.0 );
  ps_Draw_Polyline2d ( NNN-n, cp );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < NNN-n; i++ )
    ps_Mark_Circle ( (float)cp[i].x, (float)cp[i].y );

  mbs_FindBSDerivativeC2d ( n, NNN, kn, cp, NULL, dkn, dcp );
  DrawKnots ( NN-2, dkn, 1 );
  ps_Write_Command ( "100 600 translate" );
  ps_Set_Line_Width ( 6.0 );
  ps_Set_Gray ( 0.6 );
  DrawArc ( n-1, NN-2, dkn, dcp, 0.0, 2.3999 );
  DrawArc ( n-1, NN-2, dkn, dcp, 2.4, 6.0 );
  ps_Set_Line_Width ( 3.0 );
  ps_Set_Gray ( 0.0 );
  ps_Draw_Polyline2d ( NNN-n-1, dcp );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < NNN-n-1; i++ )
    ps_Mark_Circle ( (float)dcp[i].x, (float)dcp[i].y );

  ps_Fill_Circle ( 0.0, 0.0, 10.0 );
  ps_GRestore ();
} /*ShowIt*/

int main ()
{
  pkv_InitScratchMem ( 65536 );
  ps_WriteBBox ( 4, 12, 301, 302 );
  ps_OpenFile ( fn, 600 );

  ShowIt ( 0 );

  ps_CloseFile ();
  pkv_DestroyScratchMem ( );
  exit ( 0 );
} /*main*/

