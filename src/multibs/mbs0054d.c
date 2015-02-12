
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2015                            */
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

boolean mbs_multiBCDegRedd ( int ncurves, int spdimen,
                             int inpitch, int indegree, const double *inctlpoints,
                             int deltadeg,
                             int outpitch, int *outdegree, double *outctlpoints )
{
  void    *sp;
  int     n;
  double  *workspace;
  boolean result;

  if ( deltadeg > indegree )
    return false;
  else if ( !deltadeg ) {
    pkv_Selectd ( ncurves, spdimen*(indegree+1), inpitch, outpitch,
                  inctlpoints, outctlpoints );
    *outdegree = indegree;
    return true;
  }
  n = indegree-deltadeg;
  sp = workspace = pkv_GetScratchMemd ( ((n+1)*(n+2))/2 + (n+1)*(indegree+1) );
  if ( !workspace )
    return false;
  result = _mbs_multiBCDegRedd ( ncurves, spdimen, inpitch, indegree, inctlpoints,
                                 deltadeg, outpitch, outdegree, outctlpoints,
                                 workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_multiBCDegRedd*/

