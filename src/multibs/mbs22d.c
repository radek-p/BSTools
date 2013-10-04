
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

boolean mbs_BSPatchToBezd ( int spdimen,
                            int degreeu, int lastuknot, const double *uknots,
                            int degreev, int lastvknot, const double *vknots,
                            int inpitch, const double *inctlp,
                            int *kupcs, int *lastoutuknot, double *outuknots,
                            int *kvpcs, int *lastoutvknot, double *outvknots,
                            int outpitch, double *outctlp )
{
  int   NNa, MMa, pitch1, pitch2, ku, kv, skipl, skipr;
  double *ua, *va, *cpa, *cpb;
  void  *st;

  st = pkv_GetScratchMemTop ();

                /* allocate buffers */
  NNa = mbs_LastknotMaxInsd ( degreeu, lastuknot, uknots, &ku );
  MMa = mbs_LastknotMaxInsd ( degreev, lastvknot, vknots, &kv );
  pitch1 = spdimen*(lastvknot-degreev);
  pitch2 = spdimen*(MMa-degreev);
  ua = pkv_GetScratchMemd ( NNa+1 );
  va = pkv_GetScratchMemd ( MMa+1 );
  cpa = pkv_GetScratchMemd ( pitch2*(NNa-degreeu) );
  cpb = pkv_GetScratchMemd ( pitch1*(NNa-degreeu) );
  if ( !ua || !va || !cpa || !cpb )
    goto failure;

               /* pack input control points */
  pkv_Selectd ( lastuknot-degreeu, pitch1, inpitch, pitch1, inctlp, cpa );

               /* maximal "u" knot insertion */
  if ( !mbs_multiMaxKnotInsd ( 1, pitch1, degreeu, lastuknot, uknots,
                               0, cpa, &NNa, ua, 0, cpb, &skipl, &skipr ) )
    goto failure;
  NNa -= skipr+skipl;
  if ( kupcs ) *kupcs = ku;
  if ( lastoutuknot ) *lastoutuknot = NNa;
  if ( outuknots ) memmove ( outuknots, &ua[skipl], (NNa+1)*sizeof(double) );

               /* maximal "v" knot insertion */
  if ( !mbs_multiMaxKnotInsd ( NNa-degreeu, spdimen, degreev, lastvknot, vknots,
                               pitch1, &cpb[skipl*pitch1], &MMa, va, pitch2, cpa,
                               &skipl, &skipr ) )
    goto failure;
  MMa -= skipl+skipr;
  if ( kvpcs ) *kvpcs = kv;
  if ( lastoutvknot ) *lastoutvknot = MMa;
  if ( outvknots ) memmove ( outvknots, &va[skipl], (MMa+1)*sizeof(double) );

  pkv_Selectd ( NNa-degreeu, spdimen*kv*(degreev+1), pitch2, outpitch,
                &cpa[spdimen*skipl], outctlp );

  pkv_SetScratchMemTop ( st );
  return true;

failure:
  pkv_SetScratchMemTop ( st );
  return false;
} /*mbs_BSPatchToBezd*/

