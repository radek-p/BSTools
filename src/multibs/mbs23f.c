
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

void mbs_multiBSCurvesToBezf ( int spdimen, int ncurves,
                               int degree, int lastinknot, const float *inknots,
                               int inpitch, const float *inctlp,
                               int *kpcs, int *lastoutknot, float *outknots,
                               int outpitch, float *outctlp )
{
  int   NNa, auxpitch, ku, skipl, skipr;
  float *ua, *cpa;
  void  *st;

  st = pkv_GetScratchMemTop ();

                /* allocate buffers */
  NNa = mbs_LastknotMaxInsf ( degree, lastinknot, inknots, &ku );
  auxpitch = spdimen*(NNa-degree);
  ua = pkv_GetScratchMemf ( NNa+1 );
  cpa = pkv_GetScratchMemf ( auxpitch*ncurves );

               /* maximal knot insertion */
  mbs_multiMaxKnotInsf ( ncurves, spdimen, degree, lastinknot, inknots, inpitch,
                         inctlp, &NNa, ua, auxpitch, cpa, &skipl, &skipr );
  NNa -= skipr+skipl;
  if ( kpcs ) *kpcs = ku;
  if ( lastoutknot ) *lastoutknot = NNa;
  if ( outknots ) memmove ( outknots, &ua[skipl], (NNa+1)*sizeof(float) );

  pkv_Selectf ( ncurves, spdimen*(NNa-degree), auxpitch, outpitch,
                &cpa[spdimen*skipl], outctlp );

  pkv_SetScratchMemTop ( st );
} /*mbs_BSCurvesToBezf*/

