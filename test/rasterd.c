
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "psout.h"
#include "multibs.h"

char fn1[] = "rasterbcd.ps";
char fn2[] = "rasterbsd.ps";

/* ///////////////////////////////////////////////////////////////////////// */
point2d bc2[4] = {{100.0,50.0},{500.0,250.0},
                 {0.0,250.0},{400.0,50.0}};
point3d bc3[4] = {{100.0,200.0,1.0},{1000.0,420.0,2.0},
                 {0.0,450.0,1.0},{400.0,100.0,1.0}};
point2d bs2[7] = {{100.0,20.0},{250.0,50.0},{500.0,200.0},{250.0,250.0},
                  {200.0,100.0},{300.0,150.0},{400.0,20.0}};
point3d bs3[7] = {{100.0,70.0,1.0},{1000.0,400.0,4.0},{250.0,125.0,0.5},
                  {250.0,300.0,1.0},
                  {400.0,300.0,2.0},{150.0,100.0,0.5},{400.0,70.0,1.0}};
double u[11] = {0.0,0.0,0.0,0.0,0.3,0.5,0.7,1.0,1.0,1.0,1.0};

static void SetPixel ( int x, int y )
{
  char s[20];

  sprintf ( s, "%d %d _px", x, y );
  ps_Write_Command ( s );
} /*SetPixel*/

static void OutputPixels ( const xpoint *buf, int n )
{
  int i;

  for ( i = 0; i < n; i++ )
    SetPixel ( buf[i].x, buf[i].y );
} /*OutputPixels*/

int main ()
{
  pkv_InitScratchMem ( 262144 );

  ps_WriteBBox ( 59, 30, 240, 164 );
  ps_OpenFile ( fn1, 600 );
  ps_DenseScreen ();
  ps_Write_Command ( "/m 5 def" );
  ps_Write_Command ( "/_px { newpath exch m mul exch m mul moveto");
  ps_Write_Command ( "  m 0 rlineto 0 m rlineto m neg 0 rlineto closepath fill" );
  ps_Write_Command ( "} bind def" );
  ps_Set_Gray ( 0.25 );
  mbs_RasterizeBC2d ( 3, bc2, OutputPixels, true );
  ps_Set_Gray ( 0.56 );
  mbs_RasterizeBC2Rd ( 3, bc3, OutputPixels, true );
  ps_CloseFile ();

  ps_WriteBBox ( 59, 12, 240, 164 );
  ps_OpenFile ( fn2, 600 );
  ps_DenseScreen ();
  ps_Write_Command ( "/m 5 def" );
  ps_Write_Command ( "/_px { newpath exch m mul exch m mul moveto");
  ps_Write_Command ( "  m 0 rlineto 0 m rlineto m neg 0 rlineto closepath fill" );
  ps_Write_Command ( "} bind def" );
  ps_Set_Gray ( 0.25 );
  mbs_RasterizeBS2d ( 3, 10, u, bs2, OutputPixels, true );
  ps_Set_Gray ( 0.56 );
  mbs_RasterizeBS2Rd ( 3, 10, u, bs3, OutputPixels, true );
  ps_CloseFile ();

  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

