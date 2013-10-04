
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

boolean mbs_BSPatchToBezf ( int spdimen,
                            int degreeu, int lastuknot, const float *uknots,
                            int degreev, int lastvknot, const float *vknots,
                            int inpitch, const float *inctlp,
                            int *kupcs, int *lastoutuknot, float *outuknots,
                            int *kvpcs, int *lastoutvknot, float *outvknots,
                            int outpitch, float *outctlp )
{
  int   NNa, MMa, pitch1, pitch2, ku, kv, skipl, skipr;
  float *ua, *va, *cpa, *cpb;
  void  *st;

  st = pkv_GetScratchMemTop ();

                /* allocate buffers */
  NNa = mbs_LastknotMaxInsf ( degreeu, lastuknot, uknots, &ku );
  MMa = mbs_LastknotMaxInsf ( degreev, lastvknot, vknots, &kv );
  pitch1 = spdimen*(lastvknot-degreev);
  pitch2 = spdimen*(MMa-degreev);
  ua = pkv_GetScratchMemf ( NNa+1 );
  va = pkv_GetScratchMemf ( MMa+1 );
  cpa = pkv_GetScratchMemf ( pitch2*(NNa-degreeu) );
  cpb = pkv_GetScratchMemf ( pitch1*(NNa-degreeu) );
  if ( !ua || !va || !cpa || !cpb )
    goto failure;

               /* pack input control points */
  pkv_Selectf ( lastuknot-degreeu, pitch1, inpitch, pitch1, inctlp, cpa );

               /* maximal "u" knot insertion */
  if ( !mbs_multiMaxKnotInsf ( 1, pitch1, degreeu, lastuknot, uknots,
                               0, cpa, &NNa, ua, 0, cpb, &skipl, &skipr ) )
    goto failure;
  NNa -= skipr+skipl;
  if ( kupcs ) *kupcs = ku;
  if ( lastoutuknot ) *lastoutuknot = NNa;
  if ( outuknots ) memmove ( outuknots, &ua[skipl], (NNa+1)*sizeof(float) );

               /* maximal "v" knot insertion */
  if ( !mbs_multiMaxKnotInsf ( NNa-degreeu, spdimen, degreev, lastvknot, vknots,
                               pitch1, &cpb[skipl*pitch1], &MMa, va, pitch2, cpa,
                               &skipl, &skipr ) )
    goto failure;
  MMa -= skipl+skipr;
  if ( kvpcs ) *kvpcs = kv;
  if ( lastoutvknot ) *lastoutvknot = MMa;
  if ( outvknots ) memmove ( outvknots, &va[skipl], (MMa+1)*sizeof(float) );

  pkv_Selectf ( NNa-degreeu, spdimen*kv*(degreev+1), pitch2, outpitch,
                &cpa[spdimen*skipl], outctlp );

  pkv_SetScratchMemTop ( st );
  return true;

failure:
  pkv_SetScratchMemTop ( st );
  return false;
} /*mbs_BSPatchToBezf*/

