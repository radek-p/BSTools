
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "bsmesh.h"
#include "bsmprivate.h"

/* ////////////////////////////////////////////////////////////////////////// */
int bsm_HalfedgeLoopLength ( int nv, BSMvertex *mv, int *mvhei,
                             int nhe, BSMhalfedge *mhe,
                             int he )
{
  int cnt, v, vfhe, vd, e;

  if ( mhe[he].otherhalf >= 0 )
    return -1;

  cnt = 0;
  e = he;
  do {
    cnt ++;  /* count the halfedge */
    v = mhe[e].v1;
    vd = mv[v].degree;
    vfhe = mv[v].firsthalfedge;
    e = mvhei[vfhe+vd-1];
  } while ( e != he );
  return cnt;
} /*bsm_HalfedgeLoopLength*/

