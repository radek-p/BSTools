
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

void _mbs_TabBezC2Coonsd ( int spdimen, int nknu, int nknv,
                           const double *c, const double *d, const double *p,
                           const double *hu, const double *hv, double *pp,
                           double *workspace )
{
  int i, j, k, l;

  memcpy ( workspace, d, 6*nknv*spdimen*sizeof(double) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 6; j++ )
      for ( k = 0; k < 6; k++ )
        for ( l = 0; l < spdimen; l++ )
          workspace[(6*i+j)*spdimen+l] -= hv[6*i+k]*p[(6*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(double) );
  for ( i = 0; i < nknu; i++ )
    for ( j = 0; j < nknv; j++ )
      for ( k = 0; k < 6; k++ )
        for ( l = 0; l < spdimen; l++ )
          pp[(nknv*i+j)*spdimen+l] += c[(6*i+k)*spdimen+l]*hv[6*j+k] +
                                      workspace[(6*j+k)*spdimen+l]*hu[6*i+k];
} /*_mbs_TabBezC2Coonsd*/

