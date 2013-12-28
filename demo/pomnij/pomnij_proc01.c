
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <sys/types.h>
#include <sys/times.h>
#include <signal.h>
#include <unistd.h>
#include <setjmp.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <fpu_control.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "g2blendingd.h"
#include "xgedit.h"
#include "xgeipc.h"

#include "bsfile.h"

#include "pomnijipc.h"
#include "pomnij_proc.h"
#include "proc_regmem.h"

/* ////////////////////////////////////////////////////////////////////////// */
void DumpData ( void )
{
  FILE *f;
  int  i, ncp;
  point4d *cp;

  f = fopen ( "data.txt", "w+" );
        /* dump the patch */
  fprintf ( f, "  int degu = %d, degv = %d, lknu = %d, lknv = %d\n",
            degu, degv, lknu, lknv );
  fprintf ( f, "  double knotsu[] = {" );
  for ( i = 0; i < lknu; i++ )
    fprintf ( f, "%f,", knotsu[i] );
  fprintf ( f, "%f};\n", knotsu[lknu] );
  fprintf ( f, "  double knotsv[] = {" );
  for ( i = 0; i < lknv; i++ )
    fprintf ( f, "%f,", knotsv[i] );
  fprintf ( f, "%f};\n", knotsv[lknv] );
  ncp = (lknu-degu)*(lknv-degv);
  cp = (point4d*)cpoints;
  fprintf ( f, "  point3d cpoints[] = {\n" );
  for ( i = 0; i < ncp-1; i++ )
    fprintf ( f, "    {%f,%f,%f},\n", cp[i].x, cp[i].y, cp[i].z );
  fprintf ( f, "    {%f,%f,%f}};\n",  cp[ncp-1].x, cp[ncp-1].y, cp[ncp-1].z );
  fprintf ( f, "\n" );
        /* dump the constraints */
  fprintf ( f, "  int nckn = %d\n", nucknots );
  fprintf ( f, "  double ckn[] = {" );
  for ( i = 0; i < nucknots-1; i++ )
    fprintf ( f, "%f,", ucknots[i] );
  fprintf ( f, "%f};\n", ucknots[nucknots-1] );
  ncp = nucknots*(lknv-degv);
  cp = (point4d*)ucurvcp;
  fprintf ( f, "  point3d ccp[] = {\n" );
  for ( i = 0; i < ncp-1; i++ )
    fprintf ( f, "    {%f,%f,%f},\n", cp[i].x, cp[i].y, cp[i].z );
  fprintf ( f, "    {%f,%f,%f}};\n",  cp[ncp-1].x, cp[ncp-1].y, cp[ncp-1].z );
  fprintf ( f, "\n" );
        /* dump the constraint equations */
  /* ************** */
  fclose ( f );
} /*DumpData*/

void WriteThePatch ( int spdimen, int udeg, int lastknotu, double *knotsu,
                     int vdeg, int lastknotv, double *knotsv,
                     int pitch, double *cpoints )
{
  if ( bsf_OpenOutputFile ( "pp.bs", false ) ) {
    bsf_WriteBSplinePatchd ( spdimen, spdimen, false, udeg, lastknotu, knotsu,
                             vdeg, lastknotv, knotsv, false, false,
                             pitch, cpoints, NULL, NULL, NULL, NULL );
    bsf_CloseOutputFile ();
    printf ( "%s\n", "pp.bs" );
  }
} /*WriteThePatch*/

void WriteTheConstraints ( int spdimen, int nucknots, double *ucknots,
                           int vdeg, int lastknotv, double *knotsv,
                           int pitch, double *cp )
{
  int i;

  if ( nucknots <= 0 )
    return;
  if ( bsf_OpenOutputFile ( "pc.bs", false ) ) {
    for ( i = 0; i < nucknots; i++ )  
      bsf_WriteBSplineCurved ( spdimen, spdimen, false,
                               vdeg, lastknotv, knotsv, false,
                               &cp[i*pitch], NULL, NULL, NULL, NULL );
    bsf_WriteComment ( "'u' knots of interpolation" );
    bsf_WriteKnotSequenced ( nucknots-1, ucknots, false );
    bsf_CloseOutputFile ();
    printf ( "%s\n", "pc.bs" );
  }
} /*WriteTheConstraints*/

