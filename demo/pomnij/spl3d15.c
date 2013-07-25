
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkgeom.h"
#include "multibs.h"
#include "convh.h"
#include "camerad.h"
#include "xgedit.h"

#include "spl3d.h"


static void CWriteCamera ( FILE *f )
{
  fprintf ( f, "/* camera parameters:*/\n" );
  fprintf ( f, "/* position = {%f,%f,%f} */\n",
            swind.CPos[3].position.x, swind.CPos[3].position.y,
            swind.CPos[3].position.z );
  fprintf ( f, "/* psi = %f, theta = %f, phi = %f */\n",
            swind.CPos[3].psi, swind.CPos[3].theta, swind.CPos[3].phi );
  fprintf ( f, "/* f = %f */\n\n", swind.CPos[3].vd.persp.f );
} /*CWriteCamera*/

static boolean CWriteBSPatch ( FILE *f,
                               int n, int lastknotu, const double *knotsu,
                               int m, int lastknotv, const double *knotsv,
                               const point4d *cp )
{
  int i, j;

  fprintf ( f, "/*B-spline patch*/\n" );
  fprintf ( f, "int deg_u = %d;\n", n );
  fprintf ( f, "int lastknot_u = %d;\n", lastknotu );
  fprintf ( f, "double knots_u[] = {" );
  for ( i = 0, j = 3;  i <= lastknotu;  i++, j++ ) {
    fprintf ( f, "%f", knotsu[i] );
    if ( i < lastknotu )
      fprintf ( f, ",");
    else
      fprintf ( f, "};" );
    if ( j >= 7 || i == lastknotu ) {
      fprintf ( f, "\n" );
      j = 0;
    }
  }
  fprintf ( f, "int deg_v = %d;\n", m );
  fprintf ( f, "int lastknot_v = %d;\n", lastknotv );
  fprintf ( f, "double knots_v[] = {" );
  for ( i = 0, j = 3;  i <= lastknotv;  i++, j++ ) {
    fprintf ( f, "%f", knotsv[i] );
    if ( i < lastknotv )
      fprintf ( f, ",");
    else
      fprintf ( f, "};" );
    if ( j >= 7 || i == lastknotv ) {
      fprintf ( f, "\n" );
      j = 0;
    }
  }
  fprintf ( f, "point3f cp[] = {" );
  for ( i = 0, j = 2;  i < (lastknotu-n)*(lastknotv-m);  i++, j++ ) {
    fprintf ( f, "{%f,%f,%f}", cp[i].x, cp[i].y, cp[i].z );
    if ( i < (lastknotu-n)*(lastknotv-m)-1 )
      fprintf ( f, "," );
    else
      fprintf ( f, "};" );
    if ( j >= 3 || i == (lastknotu-n)*(lastknotv-m)-1 ) {
      fprintf ( f, "\n" );
      j = 0;
    }
  }
  return true;
} /*CWriteBSPatch*/

boolean CExport ( void )
{
  FILE  *f;

  f = fopen ( "pomnij.c", "w+" );
  if ( !f )
    return false;

  CWriteCamera ( f );
  fprintf ( f, "\n" );
  CWriteBSPatch ( f, degree_u, lastknot_u, knots_u,
                     degree_v, lastknot_v, knots_v, cpoints );
  fclose ( f );
  printf ( "pomnij.c\n" );
  return true;
} /*CExport*/

