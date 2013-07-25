
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"
#include "psout.h"

#define TEST7

char fn[] = "degredd.ps";

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

static void DrawGraphs ( int degree, int lastknot, const double *knots,
                         double u0, double uN,
                         double ax, double bx, double ay, double by, char dom )
{
#define DD 50
  int     i, j, fnz, nnz;
  double   t, x, y, f[10];
  point2d c[DD+1];

  ps_Set_Gray ( 0.0 );
  psl_SetLine ( (float)bx, (float)by, (float)bx, (float)(by+ay), 0.0, 1.0 );
  psl_Draw ( 0.0, 1.2, 2.0 );
  psl_Tick ( 0.0 );
  psl_Tick ( 1.0 );
  psl_Arrow ( 1.35, true );
  psl_SetLine ( (float)bx, (float)by, (float)(bx+ax), (float)by, (float)u0, (float)uN );
  psl_Draw ( (float)u0, (float)(uN+0.3), 2.0 );
  psl_Draw ( (float)knots[degree], (float)knots[lastknot-degree], 6.0 );
  psl_Arrow ( (float)(uN+0.35), true );
  t = u0-1.0e10;
  for ( i = 0; i <= lastknot; i++ ) {
    psl_Tick ( (float)knots[i] );
    if ( knots[i] > t ) {
      psl_GetPointd ( knots[i], &x, &y );
    }
    y -= 40.0;
    ps_Fill_Circle ( (float)x, (float)y, 10.0 );
    t = knots[i];
  }
  for ( i = 0; i < lastknot-degree; i++ ) {
    if ( dom ) {
      ps_Set_Gray ( 0.0 );
      ps_Write_Command ( "0 setlinecap" );
      psl_SetLine ( (float)bx, (float)(by-190.0-25*(i+1)), (float)(bx+ax),
                    (float)(by-190.0-25*(i+1)), (float)u0, (float)uN);
      psl_Draw ( (float)knots[i], (float)knots[i+degree+1], 6.0 );
      ps_Write_Command ( "1 setlinecap" );
      for ( j = i; j <= i+degree+1; j++ )
        psl_Tick ( (float)knots[j] );
    }
    ps_Set_Gray ( 0.68 );
    ps_Set_Line_Width ( 6.0 );
    for ( j = 0; j <= DD; j++ ) {
      t = knots[i] + (double)j/(double)DD*(knots[i+degree+1]-knots[i]);
      mbs_deBoorBasisd ( degree, lastknot, knots, t, &fnz, &nnz, f );
      c[j].x = bx+ax*(t-u0)/(uN-u0);
      c[j].y = by;
      if ( i >= fnz && i < fnz+nnz )
        c[j].y += ay*f[i-fnz];
    }
    ps_Draw_Polyline2d ( DD+1, c );
  }
  ps_Set_Gray ( 0.0 );
  psl_SetLine ( (float)bx, (float)by, (float)(bx+ax), (float)by, (float)u0, (float)uN );
  psl_Draw ( (float)knots[degree], (float)knots[lastknot-degree], 6.0 );
#undef DD
} /*DrawGraphs*/

#ifdef TEST1
#define deg 5
#define LKN 16
double knots[17] =
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.432665, 0.432665, 0.432665,
  0.710602, 0.710602,
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
point2d cpoints[LKN-deg] =
  {{-0.894540,-0.673759},{-0.991094,-0.299291},
   {-0.895960,0.113475},{-0.662557,0.414333},
   {-0.175038,0.651671},{0.223954,0.632677},
   {0.540992,0.532982},{0.833936,0.244086},
   {0.908158,-0.283301},{0.829225,-0.595745},
   {0.539564,-0.907801}};
#endif
#ifdef TEST2
#define deg 3
#define LKN 9
double knots[] = {
  0.0,0.0,0.0,0.0,0.455202,0.455202,1.0,1.0,1.0,1.0};
point2d cpoints[] =
{{-1.000000,0.000000},
 {-0.699173,0.232747},
 {-0.397144,0.359547},
 {0.269002,0.405357},
 {0.633640,0.278557},
 {1.000000,0.000000}};
#endif
#ifdef TEST3
#define deg  4
#define LKN 15
double knots[16] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.242069,
  0.242069, 0.435517, 0.435517, 0.735172, 0.735172, 1.000000,
  1.000000, 1.000000, 1.000000, 1.000000};
point2d cpoints[21] =
 {{-1.000000,-1.000000},
  {-1.000000,-0.757931}, {-1.000000,-0.540172},
  {-1.000000,-0.322414}, {-1.000000,-0.075862},
  {-1.000000,0.170690}, {-1.000000,0.452931},
  {-1.000000,0.735172}, {-1.000000,0.867586},
  {-1.000000,1.000000}, {-0.767921,-1.000000},
  {-0.625091,-0.885738}, {-0.547277,-0.770627},
  {-0.491743,-0.562422}, {-0.529687,-0.352997},
  {-0.596154,-0.115573}, {-0.662621,0.121850},
  {-0.738710,0.393642}, {-0.814798,0.665433},
  {-0.820927,0.812831}, {-0.767921,1.000000}};
#endif
#ifdef TEST4
#define CLOSED
#define deg  5
#define LKN 16
double knots[17] = {-0.668137,-0.343109,-0.343109,-0.191582,0.111930,
  0.111930,0.331863,0.656891,0.656891,0.808418,1.111930,1.111930,
  1.331863,1.656891,1.656891,1.808418,2.111930};
point2d cpoints[11] =
  {{0.639733,0.771523},{0.861817,0.278146},{0.570125,-0.642384},
   {-0.132587,-0.615894},{-0.997718,-0.205298},{-0.507146,0.857616},
   {0.639733,0.771523},{0.861817,0.278146},{0.570125,-0.642384},
   {-0.132587,-0.615894},{-0.997718,-0.205298}};
#endif
#ifdef TEST5
#define CLOSED
#define deg  4
#define LKN 17
double knots[18] =
  {-0.154110,-0.154110,0.000000,0.000000,0.000000,0.148973,0.148973,
   0.148973,0.845890,0.845890,0.845890,1.000000,1.000000,1.000000,
   1.148973,1.148973,1.148973,1.845890};
point2d cpoints[13] =
  {{-0.527584,-0.685418},{0.415023,-0.755457},{0.701949,-0.727013},
   {0.812485,-0.635698},{0.438558,0.239788},{0.077943,0.500792},
   {-0.335215,0.301498},{-0.903895,-0.503950},{-0.812785,-0.613058},
   {-0.527584,-0.685418},{0.415023,-0.755457},{0.701949,-0.727013},
   {0.812485,-0.635698}};
#endif
#ifdef TEST6
#define CLOSED
#define deg  4
#define LKN 18
double knots[19] =
  {-0.293948,-0.086455,-0.086455,0.118156,0.118156,0.299712,
   0.299712,0.495677,0.495677,0.706052,0.706052,0.913545,
   0.913545,1.118156,1.118156,1.299712,1.299712,1.495677,
   1.495677};
point2d cpoints[14] =
  {{-1.051984,-0.052076},{-0.888681,-0.718252},{-0.187959,-0.844698},
   {0.362558,-0.784035},{-0.159974,-0.472835},{-0.746824,-0.065027},
   {-0.191858,0.354662},{0.387734,0.718575},{-0.204909,0.801050},
   {-0.926388,0.666466},{-1.051984,-0.052076},{-0.888681,-0.718252},
   {-0.187959,-0.844698},{0.362558,-0.784035}};
#endif
#ifdef TEST7
#define CLOSED
#define deg  4
#define LKN 17
double knots[18] =
  {-0.236779,-0.236779,0.000000,0.000000,0.000000,0.384558,0.384558,
   0.384558,0.763221,0.763221,0.763221,1.000000,1.000000,1.000000,
   1.384558,1.384558,1.384558,1.763221};
point2d cpoints[13] =
  {{-0.788913,-0.491961},{0.197559,-0.590142},{0.530920,-0.502398},
   {0.587095,-0.266141},{0.148461,0.497495},{-0.108977,0.613596},
   {-0.406230,0.466930},{-0.954068,-0.198525},{-0.969273,-0.381677},
   {-0.788913,-0.491961},{0.197559,-0.590142},{0.530920,-0.502398},
   {0.587095,-0.266141}};
#endif

double outknots[50];
point2d outcpoints[50];
int outN, outdeg;

int main ( void )
{
  int i;

  pkv_InitScratchMem ( 65536 );

  for ( i = 0; i < LKN-deg; i++ ) {
    cpoints[i].x = 500.0*(cpoints[i].x+1.25);
    cpoints[i].y = 500.0*(cpoints[i].y+1.2);
  }

  ps_WriteBBox ( 34, 11, 217, 120 );
  ps_OpenFile ( fn, 600 );

  mbs_multiBSDegRedd ( 1, 2, deg, LKN, knots, 0, (double*)cpoints, 1,
                       &outdeg, &outN, outknots, 0, (double*)outcpoints );

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
#ifdef CLOSED
  ps_Write_Command ( "0 1600 translate" );
  mbs_multiBSDegRedClosedd ( 1, 2, deg, LKN, knots, 0, (double*)cpoints, 1,
                       &outdeg, &outN, outknots, 0, (double*)outcpoints );

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
#endif
  ps_CloseFile ();
  printf ( "%s\n", fn );
  printf ( "Scratch memory used: %d\n", (int)pkv_MaxScratchTaken () );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

