
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* this file was written by Mateusz Markowski                                */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "g1blendingd.h"

#include "g1blprivated.h"
#include "msgpool.h"

#define _DEBUG

/* ////////////////////////////////////////////////////////////////////////// */
/* The criterion for "lazy" recomputing of the integrals in the computations  */
/* of the Hessian matrix is still to be developed; so far the time savings    */
/* are about 1%, which is not satisfactory for me.                            */
/* /////////////////////////////////////////
czy trzeba obliczac na nowo hesjan czy nie -> czy punkty zmienily sie
///////////////////////////////// */
boolean _g1bl_LazyHessiand ( int lastknotu, int lastknotv,
                     int ni, int nsq, point3d *acp, point3d *hcp,
                     char *dirtypt, char *dirtysq, boolean *all )
{
#define THR 1.0e-14  /* threshold */
  void     *sp;
  double   *dd, diam, d;
  vector3d dcp;
  point3d  *scp;
  int      i, j, k, ip, jp, kp, nsu, nsv, pitch1;

  sp = pkv_GetScratchMemTop ();
  dd = pkv_GetScratchMemd ( ni );
  scp = pkv_GetScratchMem ( 3*3*sizeof(point3d) );
  if ( !dd || !scp ) {
    pkv_SignalError ( LIB_G1BLENDING, 13, ERRMSG_0 );
    goto failure;
  }

  pitch1 = lastknotv-DEG;
  for ( i = DEG, k = 0;  i < lastknotu-DEG2;  i++ )
    for ( j = DEG;  j < lastknotv-DEG2;  j++, k++ ) {
      kp = i*pitch1+j;
      SubtractPoints3d ( &acp[kp], &hcp[kp], &dcp );
      dd[k] = DotProduct3d ( &dcp, &dcp );
    }

  memset ( dirtysq, 0, nsq );
  nsu = lastknotu-DEG2;
  nsv = lastknotv-DEG2;
  for ( i = 0; i < nsu; i++ )
    for ( j = 0; j < nsv; j++ ) {
        /* find the diameter of the control points associated with the square */
      pkv_Selectd ( 3, 3*3, 3*pitch1, 3*3, &acp[i*pitch1+j], scp );
      diam = 0.0;
      for ( ip = 1; ip < 3*3; ip++ )
        for ( jp = 0; jp < ip; jp++ ) {
          SubtractPoints3d ( &scp[ip], &scp[jp], &dcp );
          d = DotProduct3d ( &dcp, &dcp );
          diam = max ( diam, d );
        }

      for ( ip = 0; ip < 3; ip++ )
        if ( i+ip >= DEG && i+ip < lastknotu-DEG2 )
          for ( jp = 0; jp < 3; jp++ )
            if ( j+jp >= DEG && j+jp < lastknotv-DEG2 ) {
              kp = (i+ip-DEG)*(lastknotv-DEG3) + j+jp-DEG;
              if ( dd[kp] >= diam*THR ) {
                dirtysq[i*(lastknotv-DEG2)+j] |= DIRTY_HESS;
                goto nextsq;
              }
            }
nextsq:
        ;
    }

  for ( i = k = 0; i < nsq; i++ )
    if ( dirtysq[i] )
      k ++;
  if ( (*all = k == nsq) ) {
    memcpy ( hcp, acp, (lastknotu-DEG)*(lastknotv-DEG)*sizeof(point3d) );
    memset ( dirtypt, true, ni );
    goto finish;
  }

#ifdef _DEBUG
printf ( "(%d,", k );
#endif

  memset ( dirtypt, false, ni );
  for ( ip = kp = 0;  ip < lastknotu-DEG3;  ip++ )
    for ( jp = 0;  jp < lastknotv-DEG3;  jp++, kp++ ) {
      for ( i = 0; i < 3; i++ )
        for ( j = 0; j < 3; j++ ) {
          k = (ip+i)*(lastknotv-DEG2) + jp+j;
          if ( !dirtysq[k] )
            goto nextpoint;
        }
      dirtypt[kp] = true;
nextpoint:
      ;
    }

#ifdef _DEBUG
for ( i = k = 0; i < ni; i++ )
  if ( dirtypt[i] )
    k ++;
printf ( "%d) ", k );
#endif

  for ( i = DEG, k = 0;  i < lastknotu-DEG2;  i++ )
    for ( j = DEG;  j < lastknotv-DEG2;  j++, k++ )
      if ( dirtypt[k] ) {
        kp = i*pitch1+j;
        hcp[kp] = acp[kp];
      }

finish:
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef SQ_THR
} /*_g1bl_LazyHessiand*/

