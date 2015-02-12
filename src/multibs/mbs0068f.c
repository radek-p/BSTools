
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

void _mbs_TabBezC1Coons0f ( int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp,
                   float *workspace )
{
  int i, j, k, l;

  memcpy ( workspace, d, 2*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 2; j++ )
      for ( k = 0; k < 2; k++ )
        for ( l = 0; l < spdimen; l++ )
          workspace[(2*i+j)*spdimen+l] -= hv[4*i+2*k]*p[(2*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknu; i++ )
    for ( j = 0; j < nknv; j++ )
      for ( k = 0; k < 2; k++ )
        for ( l = 0; l < spdimen; l++ )
          pp[(nknv*i+j)*spdimen+l] += c[(2*i+k)*spdimen+l]*hv[4*j+2*k] +
                                      workspace[(2*j+k)*spdimen+l]*hu[4*i+2*k];
} /*_mbs_TabBezC1Coons0f*/

