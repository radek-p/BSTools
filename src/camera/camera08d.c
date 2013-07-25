
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "camera.h"

#include "msgpool.h"


/* /////////////////////////////////////////////////////////////////////////// */
boolean CameraClipPolygon3d ( CameraRecd *CPos, int n, const point3d *p,
                              void (*output)(int n, point3d *p) )
{
  void     *sp;
  point3d  *inb, *outb, q, s;
  vector4d *cpl;
  double   tq, tp;
  int      inbufsize, outbufsize;
  int      i, j, m;

#define OUTPUTVERTEX(x) \
  { if ( m >= outbufsize ) \
      { if ( !pkv_GetScratchMem(sizeof(point3d)) ) goto failure; } \
    outbufsize++; \
    outb[m++] = x; \
  }

  sp = pkv_GetScratchMemTop ();
  if ( !(inb = pkv_GetScratchMem ( n*sizeof(point3d) )) )
    goto failure;
  inbufsize = n;
  memcpy ( inb, p, n*sizeof(point3d) );
  outb = pkv_GetScratchMemTop ();
  outbufsize = 0;

  for ( i = 0; i < CPos->ncplanes; i++ ) {
    cpl = &CPos->cplane[i];

        /* Sutherland-Hodgman clipping algorithm */
    m = 0;         /* no vertex has been output yet */
    q = inb[n-1];
    tq = DotProduct3d ( (vector3d*)cpl, &q ) + cpl->w;
    for ( j = 0; j < n; j++ ) {
      tp = DotProduct3d ( (vector3d*)cpl, &inb[j] ) + cpl->w;
      if ( tp >= 0.0 ) {
        if ( tq < 0.0 ) {
          InterPoint3d ( &q, &inb[j], -tq/(tp-tq), &s );
          OUTPUTVERTEX ( s )
        }
        OUTPUTVERTEX ( inb[j] )
      }
      else if ( tq >= 0.0 ) {
        InterPoint3d ( &q, &inb[j], -tq/(tp-tq), &s );
        OUTPUTVERTEX ( s )
      }
      q = inb[j];
      tq = tp;
    }

    if ( m > 0 ) {
      memmove ( inb, outb, m*sizeof(point3d) );
      inbufsize  += m-n;
      outbufsize += n-m;
      outb = &inb[inbufsize];
      n = m;
    }
    else
      goto success;
  }
  if ( n > 0 ) {
    for ( j = 0; j < n; j++ )
      CameraProjectPoint3d ( CPos, &inb[j], &inb[j] );
    output ( n, inb );
  }

success:
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef OUTPUTVERTEX
} /*CameraClipPolygon3d*/

