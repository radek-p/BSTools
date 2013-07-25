
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
#include "g2blendingf.h"

#define _DEBUG

/* This procedure sets up the constraint equations for a bicubic blending  */
/* patch, given by a set of curves of constant "u" parameter. The system   */
/* is in the form consistent with the optimization method assumed by the   */
/* procedures g2bl_SetupTriharmAMatrixf and g2bl_SetupTriharmRHSf, i.e. the*/
/* unknown variables are the internal control points of the patch, located */
/* in the space of dimension spdimen and the right hand side matrix has    */
/* spdimen columns.                                                        */

boolean g2bl_SetupULConstraintsf ( int lastknotu, int lastknotv, int spdimen,
                                   int ppitch, float *cp,
                                   int nucurv, float *ucknots,
                                   int cpitch, float *uccp,
                                   int *nconstr, float *cmat, float *crhs )
{
  void  *sp;
  int   nvars, ncon, bls;
  float *sknots, *w2;
  float Nfunc[4], aa;
  int   i, j, nk, fnz, nnz;

  sp = pkv_GetScratchMemTop ();
  if ( lastknotu <= 9 || lastknotv <= 9 || nucurv < 1 )
    goto failure;
        /* verify the u parameters of the constant parameter curves */
  if ( nucurv > 1 ) {
    sknots = pkv_GetScratchMemf ( nucurv );
    if ( !sknots )
      goto failure;
    memcpy ( sknots, ucknots, nucurv*sizeof(float) );
    if ( pkv_SortFast ( sizeof(float), ID_IEEE754_FLOAT, sizeof(float),
                        0, nucurv, sknots ) != SORT_OK )
      goto failure;
        /* are all knots within the range? */
    if ( sknots[0] <= 3.0 || sknots[nucurv-1] >= (float)(lastknotu-3) )
      goto failure;
        /* are all knots different? */
    for ( i = 1; i < nucurv; i++ )
      if ( sknots[i-1] >= sknots[i] )
        goto failure;
        /* verify the condition resulting from the Schoenberg-Whitney theorem */
    for ( i = 4; i < nucurv; i++ )
      if ( sknots[i-4]+4.0 >= sknots[i] )
        goto failure;
    pkv_FreeScratchMemf ( nucurv );
  }
  else {
    if ( ucknots[0] <= 3.0 || ucknots[0] >= (float)(lastknotu-3) )
      goto failure;
  }

  bls = lastknotv-9;
  nvars = (lastknotu-9)*bls;
  *nconstr = ncon = nucurv*bls;

  w2 = pkv_GetScratchMemf ( ncon*max(bls*6,nvars) );
  sknots = pkv_GetScratchMemf ( lastknotu+1 );
  if ( !w2 || !sknots )
    goto failure;
  for ( i = 0; i <= lastknotu; i++ )
    sknots[i] = (float)i;

        /* setup the constraint equations matrix */
  memset ( cmat, 0, nvars*ncon*sizeof(float) );
  memset ( w2, 0, 6*bls*ncon*sizeof(float) );
  for ( nk = 0; nk < nucurv; nk++ ) {
    mbs_deBoorBasisf ( 3, lastknotu, sknots, ucknots[nk], &fnz, &nnz, Nfunc );
    for ( i = fnz; i < 3 && i < fnz+nnz; i++ ) {
      aa = Nfunc[i-fnz];
      for ( j = 0; j < bls; j++ )
        w2[(nk*bls+j)*6*bls + i*bls + j] = aa;
    }
    for ( i = max(fnz,lastknotu-6); i < lastknotu-3 && i < fnz+nnz; i++ ) {
      aa = Nfunc[i-fnz];
      for ( j = 0; j < bls; j++ )
        w2[(nk*bls+j)*6*bls + (i-lastknotu+9)*bls + j] = aa;
    }
    for ( i = max(fnz,3); i < lastknotu-6 && i < fnz+nnz; i++ ) {
      aa = Nfunc[i-fnz];
      for ( j = 0; j < bls; j++ )
        cmat[(nk*bls+j)*nvars + (i-3)*bls + j] = aa;
    }
  }
        /* setup the constraint equations right hand side */
  pkv_Selectf ( nucurv, spdimen*bls, cpitch, spdimen*bls,
                &uccp[3*spdimen], crhs );
  for ( i = 0; i < 3; i++ )
    pkn_MultMatrixSubf ( ncon, bls, 6*bls, &w2[i*bls],
                         spdimen, spdimen, &cp[i*ppitch+3*spdimen],
                         spdimen, crhs );
  for ( i = 0; i < 3; i++ )
    pkn_MultMatrixSubf ( ncon, bls, 6*bls, &w2[(i+3)*bls],
                         spdimen, spdimen, &cp[(lastknotu-6+i)*ppitch+3*spdimen],
                         spdimen, crhs );
  
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2bl_SetupULConstraintsf*/

