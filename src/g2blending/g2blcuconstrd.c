
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

/* This procedure sets up the constraint equations for a closed bicubic    */
/* blending patch, given by a set of curves of constant "u" parameter.     */
/* The system is in the form consistent with the optimization method       */
/* assumed by the procedures g2bl_SetupClosedTriharmAMatrix and            */
/* g2bl_SetupClosedTriharmRHSd, i.e. the unknown variables are the         */
/* internal control points of the patch, located in the space of dimension */
/* spdimen and the right hand side matrix has spdimen columns.             */

/* ******************** */
/* Niedokonczone !!!!!! */
/* ******************** */
boolean g2bl_SetupClosedULConstraintsd ( int lastknotu, int lastknotv, int spdimen,
                                         int ppitch, double *cp,
                                         int nucurv, double *ucknots,
                                         int cpitch, double *uccp,
                                         int *nconstr, double *cmat, double *crhs )
{
  void   *sp;
  int    nvars, ncon, bls;
  double *sknots;
  double Nfunc[4], K;
  int    i, ii, j, nk, fnz, nnz;

  sp = pkv_GetScratchMemTop ();
  if ( lastknotu <= 9 || lastknotv <= 9 || nucurv < 1 )
    goto failure;
        /* verify the u parameters of the constant parameter curves */
  if ( nucurv > 1 ) {
    K = (double)(lastknotu-6);
    sknots = pkv_GetScratchMemd ( nucurv+3 );
    if ( !sknots )
      goto failure;
    memcpy ( sknots, ucknots, nucurv*sizeof(double) );
    if ( pkv_SortFast ( sizeof(double), ID_IEEE754_DOUBLE, sizeof(double),
                        0, nucurv, sknots ) != SORT_OK )
      goto failure;
        /* are all knots within the range? */
    if ( sknots[0] <= 3.0 || sknots[nucurv-1] >= (double)(lastknotu-3) )
      goto failure;
        /* are all knots different? */
    for ( i = 1; i < nucurv; i++ )
      if ( sknots[i-1] >= sknots[i] )
        goto failure;
        /* verify the condition resulting from the Schoenberg-Whitney theorem */
    for ( i = 0; i < 3; i++ )
      sknots[i+nucurv] = sknots[i] + K;
    for ( i = 4; i < nucurv+3; i++ )
      if ( sknots[i-4]+4.0 >= sknots[i] )
        goto failure;
    pkv_FreeScratchMemd ( nucurv+3 );
  }
  else {
    if ( ucknots[0] <= 3.0 || ucknots[0] >= (double)(lastknotu-3) )
      goto failure;
  }

  bls = lastknotv-9;
  nvars = (lastknotu-6)*bls;
  *nconstr = ncon = nucurv*bls;

  sknots = pkv_GetScratchMemd ( lastknotu+1 );
  if ( !sknots )
    goto failure;
  for ( i = 0; i <= lastknotu; i++ )
    sknots[i] = (double)i;

        /* setup the constraint equations matrix */
  memset ( cmat, 0, nvars*ncon*sizeof(double) );
  for ( nk = 0; nk < nucurv; nk++ ) {
    mbs_deBoorBasisd ( 3, lastknotu, sknots, ucknots[nk], &fnz, &nnz, Nfunc );
    for ( i = fnz; i < fnz+nnz; i++ ) {
      if ( i < lastknotu-6 )
        ii = i;
      else
        ii = i - lastknotu+6;
      for ( j = 0; j < bls; j++ )
        cmat[(nk*bls+j)*nvars + ii*bls + j] = Nfunc[i-fnz];
    }
    
  }
        /* setup the constraint equations right hand side */
  pkv_Selectd ( nucurv, spdimen*bls, cpitch, spdimen*bls,
                &uccp[3*spdimen], crhs );
  
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2bl_SetupClosedULConstraintsd*/

