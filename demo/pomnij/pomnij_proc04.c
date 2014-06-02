
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

#include "bsfile.h"

#include "pomnijipc.h"
#include "pomnij_proc.h"
#include "proc_regmem.h"

/* ////////////////////////////////////////////////////////////////////////// */
void SetupG1BLKnots ( void )
{
  int i;

  clu = knotsu[1] >= knotsu[2] && knotsu[lknu-2] >= knotsu[lknu-1];
  mbs_TransformAffKnotsd ( 2, lknu, knotsu, (double)2, (double)(lknu-2), knotsu );
  for ( i = 0; i <= lknu; i++ )
    knotsu[i] = (double)i;
  if ( clu ) {
    knotsu[0] = knotsu[1] = knotsu[2];
    knotsu[lknu] = knotsu[lknu-1] = knotsu[lknu-2];
  }

  clv = knotsv[1] >= knotsv[2] && knotsv[lknv-2] >= knotsv[lknv-1];
  mbs_TransformAffKnotsd ( 2, lknv, knotsv, (double)2, (double)(lknv-2), knotsv );
  for ( i = 0; i <= lknv; i++ )
    knotsv[i] = (double)i;
  if ( clv ) {
    knotsv[0] = knotsv[1] = knotsv[2];
    knotsv[lknv] = knotsv[lknv-1] = knotsv[lknv-2];
  }
} /*SetupG1BLKnots*/

void G1FreeBoundaryToClamped ( point3d *ccp )
{
  void   *sp;
  double newkn[3];
  int    i, pitch;

  sp = pkv_GetScratchMemTop ();
  pitch = 3*(lknv-degv);
  newkn[0] = 0.0;
  for ( i = 1; i <= 2; i++ )
    newkn[i] = 2.0;
  mbs_multiBSChangeLeftKnotsd ( lknu-degu, 3, degv, knotsv,
                                pitch, (double*)ccp, newkn );
  for ( i = 0; i < 2; i++ )
    newkn[i] = (double)(lknv-2);
  newkn[2] = (double)lknv;
  mbs_multiBSChangeRightKnotsd ( lknu-degu, 3, degv, lknv, knotsv,
                                 pitch, (double*)ccp, newkn );
  newkn[0] = 0.0;
  for ( i = 1; i <= 2; i++ )
    newkn[i] = 2.0;
  mbs_multiBSChangeLeftKnotsd ( 1, pitch, degu,
                                knotsu, 0, (double*)ccp, newkn );
  for ( i = 0; i < 2; i++ )
    newkn[i] = (double)(lknu-2);
  newkn[2] = (double)lknu;
  mbs_multiBSChangeRightKnotsd ( 1, pitch, degu, lknu,
                                 knotsu, 0, (double*)ccp, newkn );

  pkv_SetScratchMemTop ( sp );
} /*G1FreeBoundaryToClamped*/

void G1ClampedBoundaryToFree ( point3d *ccp )
{
  void   *sp;
  double newkn[3];
  int    i, pitch;

  sp = pkv_GetScratchMemTop ();
  pitch = 3*(lknv-degv);
  for ( i = 0; i <= 2; i++ )
    newkn[i] = (double)i;
  mbs_multiBSChangeLeftKnotsd ( 1, pitch, degu,
                                knotsu, 0, (double*)ccp, newkn );
  for ( i = 0; i <= 2; i++ )
    newkn[i] = (double)(lknu-2+i);
  mbs_multiBSChangeRightKnotsd ( 1, pitch, degu, lknu,
                                 knotsu, 0, (double*)ccp, newkn );
  for ( i = 0; i <= 2; i++ )
    newkn[i] = (double)i;
  mbs_multiBSChangeLeftKnotsd ( lknu-degu, 3, degv, knotsv,
                                pitch, (double*)ccp, newkn );
  for ( i = 0; i <= 2; i++ )
    newkn[i] = (double)(lknv-2+i);
  mbs_multiBSChangeRightKnotsd ( lknu-degu, 3, degv, lknv, knotsv,
                                 pitch, (double*)ccp, newkn );

  pkv_SetScratchMemTop ( sp );
} /*G1ClampedBoundaryToFree*/

boolean G1AdjustOptRange ( void )
{

printf ( "%d, %d, %d, %d\n", options.opt_range[0], options.opt_range[1],
             options.opt_range[2], options.opt_range[3] );

  options.opt_range[0] = max ( options.opt_range[0], 2 );
  options.opt_range[1] = min ( options.opt_range[1], lknu-5 );
  options.opt_range[2] = max ( options.opt_range[2], 2 );
  options.opt_range[3] = min ( options.opt_range[3], lknv-5 );
  return options.opt_range[0] <= options.opt_range[1] &&
         options.opt_range[2] <= options.opt_range[3];
} /*G1AdjustOptRange*/

/* ////////////////////////////////////////////////////////////////////////// */
void SetupG2BLKnots ( void )
{
  int i;

  clu = knotsu[1] >= knotsu[2] && knotsu[2] >= knotsu[3] &&
         knotsu[lknu-3] >= knotsu[lknu-2] && knotsu[lknu-2] >= knotsu[lknu-1];
  mbs_TransformAffKnotsd ( 3, lknu, knotsu, (double)3, (double)(lknu-3), knotsu );
  for ( i = 0; i <= lknu; i++ )
    knotsu[i] = (double)i;
  if ( clu ) {
    knotsu[0] = knotsu[1] = knotsu[2] = knotsu[3];
    knotsu[lknu] = knotsu[lknu-1] = knotsu[lknu-2] = knotsu[lknu-3];
  }

  clv = knotsv[1] >= knotsv[2] && knotsv[2] >= knotsv[3] &&
         knotsv[lknv-3] >= knotsv[lknv-2] && knotsv[lknv-2] >= knotsv[lknv-1];
  mbs_TransformAffKnotsd ( 3, lknv, knotsv, (double)3, (double)(lknv-3), knotsv );
  for ( i = 0; i <= lknv; i++ )
    knotsv[i] = (double)i;
  if ( clv ) {
    knotsv[0] = knotsv[1] = knotsv[2] = knotsv[3];
    knotsv[lknv] = knotsv[lknv-1] = knotsv[lknv-2] = knotsv[lknv-3];
  }
} /*SetupG2BLKnots*/

void G2FreeBoundaryToClamped ( point3d *ccp )
{
  void   *sp;
  double newkn[4];
  int    i, pitch;

  sp = pkv_GetScratchMemTop ();
  pitch = 3*(lknv-degv);
  newkn[0] = 0.0;
  for ( i = 1; i <= 3; i++ )
    newkn[i] = 3.0;
  mbs_multiBSChangeLeftKnotsd ( lknu-degu, 3, degv, knotsv,
                                pitch, (double*)ccp, newkn );
  for ( i = 0; i < 3; i++ )
    newkn[i] = (double)(lknv-3);
  newkn[3] = (double)lknv;
  mbs_multiBSChangeRightKnotsd ( lknu-degu, 3, degv, lknv, knotsv,
                                 pitch, (double*)ccp, newkn );
  newkn[0] = 0.0;
  for ( i = 1; i <= 3; i++ )
    newkn[i] = 3.0;
  mbs_multiBSChangeLeftKnotsd ( 1, pitch, degu,
                                knotsu, 0, (double*)ccp, newkn );
  for ( i = 0; i < 3; i++ )
    newkn[i] = (double)(lknu-3);
  newkn[3] = (double)lknu;
  mbs_multiBSChangeRightKnotsd ( 1, pitch, degu, lknu,
                                 knotsu, 0, (double*)ccp, newkn );

  pkv_SetScratchMemTop ( sp );
} /*G2FreeBoundaryToClamped*/

void G2ClampedBoundaryToFree ( point3d *ccp )
{
  void   *sp;
  double newkn[4];
  int    i, pitch;

  sp = pkv_GetScratchMemTop ();
  pitch = 3*(lknv-degv);
  for ( i = 0; i <= 3; i++ )
    newkn[i] = (double)i;
  mbs_multiBSChangeLeftKnotsd ( 1, pitch, degu,
                                knotsu, 0, (double*)ccp, newkn );
  for ( i = 0; i <= 3; i++ )
    newkn[i] = (double)(lknu-3+i);
  mbs_multiBSChangeRightKnotsd ( 1, pitch, degu, lknu,
                                 knotsu, 0, (double*)ccp, newkn );
  for ( i = 0; i <= 3; i++ )
    newkn[i] = (double)i;
  mbs_multiBSChangeLeftKnotsd ( lknu-degu, 3, degv, knotsv,
                                pitch, (double*)ccp, newkn );
  for ( i = 0; i <= 3; i++ )
    newkn[i] = (double)(lknv-3+i);
  mbs_multiBSChangeRightKnotsd ( lknu-degu, 3, degv, lknv, knotsv,
                                 pitch, (double*)ccp, newkn );

  pkv_SetScratchMemTop ( sp );
} /*G2ClampedBoundaryToFree*/

boolean G2AdjustOptRange ( void )
{
  options.opt_range[0] = max ( options.opt_range[0], 3 );
  options.opt_range[1] = min ( options.opt_range[1], lknu-7 );
  options.opt_range[2] = max ( options.opt_range[2], 3 );
  options.opt_range[3] = min ( options.opt_range[3], lknv-7 );
  return options.opt_range[0] <= options.opt_range[1] &&
         options.opt_range[2] <= options.opt_range[3];
} /*G2AdjustOptRange*/

/* ////////////////////////////////////////////////////////////////////////// */
void PreTransformation ( int npcp, point3d *pcp )
{
  int i;

  for ( i = 0; i < npcp; i++ )
    TransPoint3d ( &options.trans, &pcp[i], &pcp[i] );
} /*PreTransformation*/

boolean PostTransformation ( int npcp, point3d *pcp )
{
  trans3d itr;
  int     i;

  itr = options.trans;
  if ( InvertTrans3d ( &itr ) ) {
    for ( i = 0; i < npcp; i++ )
      TransPoint3d ( &itr, &pcp[i], &pcp[i] );
    return true;
  }
  else
    return false;
} /*PostTransformation*/

