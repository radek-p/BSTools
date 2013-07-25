
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

boolean mbs_multiLaneRiesenfeldd ( int spdimen, int ncurves, int degree,
                         int inlastknot, int inpitch, const double *incp,
                         int *outlastknot, int outpitch, double *outcp )
{
  void   *sp;
  int    nn, i, j, k, l, p, q;
  double *acp;

  sp = pkv_GetScratchMemTop ();
  if ( ncurves < 1 || spdimen < 1 || degree < 1 || inlastknot < 2*degree )
    goto failure;
  nn = (2*(inlastknot-degree)-1);
  acp = pkv_GetScratchMemd ( nn*spdimen );
  if ( !acp )
    goto failure;
  *outlastknot = 2*(inlastknot-degree);

        /* Lane-Riesenfeld algorithm */
  for ( k = p = q = 0;  k < ncurves;  k++, p += inpitch, q += outpitch ) {
          /* refining */
    for ( i = 0; i < inlastknot-degree-1; i++ ) {
      memcpy ( &acp[2*i*spdimen], &incp[p+i*spdimen], spdimen*sizeof(double) );
      for ( l = 0; l < spdimen; l++ )
        acp[(2*i+1)*spdimen+l] =
                0.5*(incp[p+i*spdimen+l]+incp[p+(i+1)*spdimen+l]);
    }
    memcpy ( &acp[2*i*spdimen], &incp[p+i*spdimen], spdimen*sizeof(double) );
          /* averaging */
    for ( j = 1; j < degree; j++ )
      for ( i = 0; i < 2*(inlastknot-degree)-1-j; i++ )
        for ( l = 0; l < spdimen; l++ )
          acp[i*spdimen+l] = 0.5*(acp[i*spdimen+l]+acp[(i+1)*spdimen+l]);
          /* storing the result at the destination place */
    memcpy ( &outcp[q], acp, ((*outlastknot)-degree)*spdimen*sizeof(double) );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiLaneRiesenfeldd*/

