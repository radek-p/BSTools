
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "g2blendingd.h"

#include "g2blprivated.h"

#define _DEBUG

/* ////////////////////////////////////////////////////////////////////////// */
/* The criterion for "lazy" recomputing of the integrals in the computations  */
/* of the Hessian matrix is still to be developed; so far the time savings    */
/* are about 1%, which is not satisfactory for me.                            */
/* ////////////////////////////////////////////////////////////////////////// */
boolean _g2bl_ClosedLazyHessiand ( int lastknotu, int lastknotv,
                     int niip, int nsq, point3d *acp, point3d *hcp,
                     char *dirtypt, char *dirtysq, boolean *all )
{
#define THR 1.0e-14  /* threshold */
  void     *sp;
  double   *dd, diam, d;
  vector3d dcp;
  point3d  *scp;
  int      i, j, k, ip, jp, kp, nsu, nsv, pitch1;

  sp = pkv_GetScratchMemTop ();
  dd = pkv_GetScratchMemd ( niip );
  scp = pkv_GetScratchMem ( 16*sizeof(point3d) );
  if ( !dd || !scp ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }

  pitch1 = lastknotv-3;
  for ( i = 0, k = 0;  i < lastknotu-3;  i++ )
    for ( j = 3;  j < lastknotv-6;  j++, k++ ) {
      kp = i*pitch1+j;
      SubtractPoints3d ( &acp[kp], &hcp[kp], &dcp );
      dd[k] = DotProduct3d ( &dcp, &dcp );
    }

  memset ( dirtysq, 0, nsq );
  nsu = lastknotu-6;
  nsv = lastknotv-6;
  for ( i = 0; i < nsu; i++ )
    for ( j = 0; j < nsv; j++ ) {
        /* find the diameter of the control points associated with the square */
      pkv_Selectd ( 4, 3*4, 3*pitch1, 3*4, &acp[i*pitch1+j], scp );
      diam = 0.0;
      for ( ip = 1; ip < 16; ip++ )
        for ( jp = 0; jp < ip; jp++ ) {
          SubtractPoints3d ( &scp[ip], &scp[jp], &dcp );
          d = DotProduct3d ( &dcp, &dcp );
          diam = max ( diam, d );
        }

      for ( ip = 0; ip < 4; ip++ )
        if ( i+ip >= 0 && i+ip < lastknotu-6 )
          for ( jp = 0; jp < 4; jp++ )
            if ( j+jp >= 3 && j+jp < lastknotv-6 ) {
              kp = (i+ip)*(lastknotv-9) + j+jp-3;
              if ( dd[kp] >= diam*THR ) {
                dirtysq[i*(lastknotv-6)+j] |= DIRTY_HESS;
                goto nextsq;
              }
            }
nextsq:
        ;
    }

  for ( i = k = 0; i < nsq; i++ )
    if ( dirtysq[i] )
      k ++;

/* do poprawienia bledu nizej */
k = nsq;

  if ( (*all = k == nsq) ) {
    memcpy ( hcp, acp, (lastknotu-3)*(lastknotv-3)*sizeof(point3d) );
    memset ( dirtypt, true, niip );
    goto finish;
  }

/* nizej jest blad do poprawienia */
#ifdef _DEBUG
printf ( "(%d,", k );
#endif

  memset ( dirtypt, false, niip );
  for ( ip = kp = 0;  ip < lastknotu-3;  ip++ )
    for ( jp = 0;  jp < lastknotv-9;  jp++, kp++ ) {
      for ( i = 0; i < 4; i++ )
        for ( j = 0; j < 4; j++ ) {
          k = (ip+i)*(lastknotv-6) + jp+j;
          if ( !dirtysq[k] )
            goto nextpoint;
        }
      dirtypt[kp] = true;
nextpoint:
      ;
    }

#ifdef _DEBUG
for ( i = k = 0; i < niip; i++ )
  if ( dirtypt[i] )
    k ++;
printf ( "%d) ", k );
#endif

  for ( i = 0, k = 0;  i < lastknotu-3;  i++ )
    for ( j = 3;  j < lastknotv-6;  j++, k++ )
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
} /*_g2bl_ClosedLazyHessiand*/

