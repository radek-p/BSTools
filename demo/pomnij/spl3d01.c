
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

/* ///////////////////////////////////////////////////////////////////////// */
void ErrorHandler ( int module, const char *file, int line,
                    int errcode, const char *errstr )
{
  fprintf ( stderr, "Error in module %d, file %s, line %d: %s\n",
                     module, file, line, errstr );
  DumpData ();
  exit ( 1 );
} /*ErrorHandler*/

/* ////////////////////////////////////////////////////////////////////////// */
void DumpData ( void )
{
  FILE *f;
  int i, np;

  f = fopen ( "pomnij.dat", "w+" );
  fprintf ( f, "int degree_u = %d, lastknot_u=%d;\n",
            degree_u, lastknot_u );
  fprintf ( f, "doub;e knots_u[%d] =\n", lastknot_u+1 );
  for ( i = 0; i <= lastknot_u; i++ ) {
    if ( i == 0 ) fprintf ( f, " {" );
    fprintf ( f, "%f", knots_u[i] );
    if ( i < lastknot_u ) {
      if ( i % 6 != 5 ) fprintf ( f, ", " );
      else fprintf ( f, ",\n  ");
    }
    else fprintf ( f, "};\n" );
  }
  fprintf ( f, "int degree_v = %d, lastknot_v =  %d;\n",
            degree_v, lastknot_v );
  fprintf ( f, "double knots_v[%d] =\n", lastknot_v+1 );
  for ( i = 0; i <= lastknot_v; i++ ) {
    if ( i == 0 ) fprintf ( f, " {" );
    fprintf ( f, "%f", knots_v[i] );
    if ( i < lastknot_v ) {
      if ( i % 6 != 5 ) fprintf ( f, ", " );
      else fprintf ( f, ",\n  ");
    }
    else fprintf ( f, "};\n" );
  }
  np = (lastknot_u-degree_u)*(lastknot_v-degree_v);
  fprintf ( f, "point4f cpoints[%d] =\n", np );
  for ( i = 0; i < np; i++ ) {
    if ( i == 0 ) fprintf ( f, " {");
    fprintf ( f, "{%f,%f,%f,%f}",
              cpoints[i].x, cpoints[i].y, cpoints[i].z, cpoints[i].w );
    if ( i < np-1 ) {
      if ( i % 2 != 1 ) fprintf ( f, ", " );
      else fprintf ( f, ",\n  " );
    }
    else fprintf ( f, "};\n" );
  }
  fprintf ( f, "/*closed_u=%d, closed_v=%d*/\n",
            kwind.closed_u, kwind.closed_v );
  fclose ( f );
} /*DumpData*/

