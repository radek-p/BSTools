
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "psout.h"

char fn[] = "clcuintf.ps";

#define LIKN 4
point2f x[LIKN+1] = {{500,300},{900,300},{1050,850},{400,700},{500,300}};
float   ikn[LIKN+1] = {0.0, 1.0, 3.0, 4.0, 5.5};

point2f cp[100];
float   kn[100];
int     lastknot;

static void DrawCPolygon ( int degree, int lastknot, const point2f *d )
{
  int i;

  ps_Draw_Polyline2f ( lastknot-degree, d );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < lastknot-degree; i++ )
    ps_Mark_Circle ( d[i].x, d[i].y );
} /*DrawCPolygon*/

static void DrawBSCurve ( int degree, int lastknot, const float *u,
                          const point2f *d )
{
#define DD 200
  int     i;
  float   t;
  point2f *c;

  c = pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  for ( i = 0; i <= DD; i++ ) {
    t = u[degree] + (float)i/(float)DD*(u[lastknot-degree]-u[degree]);
    mbs_deBoorC2f ( degree, lastknot, u, d, t, &c[i] );
  }
  ps_Draw_Polyline2f ( DD+1, c );
  pkv_FreeScratchMem ( (DD+1)*sizeof(point2f) );
#undef DD
} /*DrawBSCurve*/

int main ( void )
{
  int i;

  pkv_InitScratchMem ( 1048576 );
  ps_WriteBBox ( 0, 0, 200, 180 );
  ps_OpenFile ( fn, 600 );

  mbs_multiBSCubicClosedInterpf ( LIKN, ikn, 1, 2, 0, (float*)x,
                                  &lastknot, kn, 0, (float*)cp );
  ps_Set_RGB ( 1.0, 0.0, 0.0 );
  DrawCPolygon ( 3, lastknot, cp );
  ps_Set_RGB ( 0.0, 1.0, 0.0 );
  DrawBSCurve ( 3, lastknot, kn, cp );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i <= LIKN; i++ )
    ps_Fill_Circle ( x[i].x, x[i].y, 10.0 );
  ps_CloseFile ();
  printf ( "%s\n", fn );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

