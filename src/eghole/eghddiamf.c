
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2009                            */
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

#undef CONST_
#define CONST_

#include "eg1holef.h"
#include "eghprivatef.h"
#include "eg1herror.h"

/* ///////////////////////////////////////////////////////////////////////// */
/* compute the diameter of domain, or rather the diameter of the set of */
/* Bezier control points of the domain boundary curves */
float gh_DomainDiamf ( GHoleDomainf *domain )
{
  void     *sp;
  GHolePrivateRecf *privateG;
  int      hole_k, i, j;
  point2f  *p, *q;
  vector2f v;
  float    diam, d;

  sp = pkv_GetScratchMemTop ();
  if ( !(privateG = domain->privateG) )
    goto failure;
  hole_k = domain->hole_k;
  p = (point2f*)pkv_GetScratchMem ( 6*hole_k*sizeof(point2f) );
  if ( !p )
    goto failure;
  q = privateG->surrpc;
  for ( i = j = 0;  i < hole_k;  i++ ) {
    p[j++] = q[24*i];     p[j++] = q[24*i+1];   p[j++] = q[24*i+2];
    p[j++] = q[24*i+13];  p[j++] = q[24*i+14];  p[j++] = q[24*i+15];
  }
  diam = 0.0;
  for ( i = 1; i < 6*hole_k; i++ )
    for ( j = 0; j < i; j++ ) {
      SubtractPoints2f ( &p[i], &p[j], &v );
      d = v.x*v.x+v.y*v.y;
      diam = max ( diam, d );
    }
  privateG->diam = diam = (float)sqrt ( diam );
  pkv_SetScratchMemTop ( sp );
  return diam;

failure:
  pkv_SetScratchMemTop ( sp );
  return 0.0;
} /*gh_DomainDiamf*/

