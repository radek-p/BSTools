
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

/* ////////////////////////////////////////// */
/* conversion of a B-spline patch to piecewise Bezier form */

boolean mbs_multiBSCurvesToBezd ( int spdimen, int ncurves,
                                  int degree, int lastinknot, const double *inknots,
                                  int inpitch, const double *inctlp,
                                  int *kpcs, int *lastoutknot, double *outknots,
                                  int outpitch, double *outctlp )
{
  int   NNa, auxpitch, ku, skipl, skipr;
  double *ua, *cpa;
  void  *st;

  st = pkv_GetScratchMemTop ();

                /* allocate buffers */
  NNa = mbs_LastknotMaxInsd ( degree, lastinknot, inknots, &ku );
  auxpitch = spdimen*(NNa-degree);
  ua = pkv_GetScratchMemd ( NNa+1 );
  cpa = pkv_GetScratchMemd ( auxpitch*ncurves );
  if ( !ua || !cpa )
    goto failure;

               /* maximal knot insertion */
  if ( !mbs_multiMaxKnotInsd ( ncurves, spdimen, degree,
                               lastinknot, inknots, inpitch,
                               inctlp, &NNa, ua, auxpitch, cpa,
                               &skipl, &skipr ) )
    goto failure;
  NNa -= skipr+skipl;
  if ( kpcs ) *kpcs = ku;
  if ( lastoutknot ) *lastoutknot = NNa;
  if ( outknots ) memmove ( outknots, &ua[skipl], (NNa+1)*sizeof(double) );

  pkv_Selectd ( ncurves, spdimen*(NNa-degree), auxpitch, outpitch,
                &cpa[spdimen*skipl], outctlp );

  pkv_SetScratchMemTop ( st );
  return true;

failure:
  pkv_SetScratchMemTop ( st );
  return false;
} /*mbs_BSCurvesToBezd*/

