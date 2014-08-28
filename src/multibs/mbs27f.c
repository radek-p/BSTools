
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


/* /////////////////////////////////////////// */
/* Trimmed patch boundary compilation          */

/* ///////////////////////////////////////////////////////////////////////// */
/* Trimmed patch domain boundary is converted into a code to be processed    */
/* by the procedures drawing the patch. Each command consists of one byte    */
/* opcode followed by data. The opcode may have the values                   */
/* 0,1-polyline. Followed by one byte space dimension (either 2 or 3)        */
/*     and short number N of the following line segments (N+1 vertices, i.e. */
/*     point2f or vector3f structures, depending on dimension.               */
/*     If the dimension is 2 then these specify cartesian coordinates,       */
/*     and if 3 then homogeneous coordinates of each point.                  */
/* 2,3-Bezier curve. Followed by one byte space dimension and short number n */
/*     indicating the degree. This is followed by n+1 control points, i.e.   */
/*     point2f or vector3f structures as above.                              */
/* 4 - end of code. This opcode is not followed by any further data.         */
/* The opcodes 0 and 2 indicate that the last point of the polyline or curve */
/* is linked to the first point of the next polyline curve to form a closed  */
/* path. The opcodes 1 and 3 indicate a closing element, whose end is linked */
/* to the first point of the first element after the last closing one or the */
/* first one at all.                                                         */
/* ///////////////////////////////////////////////////////////////////////// */

int mbs_TrimCVBoundSizef ( int nelem, const mbs_polycurvef *bound )
{
  int   i;
  int   size;
  int   dim, deg, N, kn;
  float *knots;

/* compute the total size of data and opcodes */

  size = 0;
  for ( i = 0; i < nelem; i++ ) {
    dim = bound[i].spdimen;
    deg = bound[i].degree;
    if ( dim < 2 || dim > 3 || deg < 1 )
      return 0;

    N   = bound[i].lastknot;
                                   /* recognize the case */
    if ( (knots = bound[i].knots) ) {  /* a B-spline curve */
      kn = mbs_NumKnotIntervalsf ( deg, N, knots );
      size += kn*(dim*(deg+1)*sizeof(float) + 2 + sizeof(short));
    }
    else if ( N == -1 )              /* a Bezier curve */
      size += dim*(deg+1)*sizeof(float) + 2 + sizeof(short);
    else if ( deg == 1 && N > 0 )    /* a polyline */
      size += dim*(N+1)*sizeof(float) + 2 + sizeof(short);
    else
      return 0;
  }
  return size + 1;
} /*mbs_TrimCVBoundSizef*/

void *mbs_CompileTrimPatchBoundf ( int nelem, const mbs_polycurvef *bound,
                                   void *buffer )
{
  void  *sp;
  int   size = 0, inpitch, outpitch;
  char  *bufp, *bp;
  short *bn;
  int   i, j;
  int   dim, deg, N, kpcs;
  float *knots, *points, *bpoints;;

  sp = pkv_GetScratchMemTop ();
  if ( !buffer ) {
    size = mbs_TrimCVBoundSizef ( nelem, bound );
    if ( !size )
      goto failure;
    buffer = pkv_GetScratchMem ( size );
    if ( !buffer )
      goto failure;
  }

  bufp = bp = (char*)buffer;
  for ( i = 0; i < nelem; i++ ) {
    bn      = (short*)(bufp+2);
    bpoints = (float*)(bufp+2+sizeof(short));

    bufp[1] = (char)(dim = bound[i].spdimen);
    deg = bound[i].degree;
    if ( dim < 2 || dim > 3 || deg < 1 )
      goto failure;
    N       = bound[i].lastknot;
    points  = bound[i].points;
                                   /* recognize the case */
    if ( (knots = bound[i].knots) ) {  /* a B-spline curve */
                                       /* convert it into a sequence */
                                       /* of Bezier curves */
      if ( !mbs_multiBSCurvesToBezf ( dim, 1, deg, N, knots, 0, points,
                                      &kpcs, NULL, NULL, 0, bpoints ) )
        goto failure;

            /* rearrange data - cannot use pkv_Rearrangef, as this is not */
            /* an array of floats, though it might work being tricky      */
      inpitch  = dim*(deg+1)*sizeof(float);
      outpitch = inpitch+2+sizeof(short);
      pkv_Rearrangec ( kpcs, inpitch, inpitch, outpitch, (char*)bpoints );

            /* insert the Bezier curve opcodes */
      for ( j = 0;  j < kpcs;  j++ ) {
        bufp[0] = 2;
        bufp[1] = (char)dim;
        bn  = (short*)(&bufp[2]);
        *bn = (short)deg;
        bp = bufp;
        bufp += outpitch;
      }
      if ( bound[i].closing )
        bp[0] = 3;
    }
    else if ( N == -1 ) {            /* a Bezier curve */
      if ( bound[i].closing )
        bufp[0] = 3;
      else
        bufp[0] = 2;
      bufp[1] = (char)dim;
      *bn = (short)deg;
      memcpy ( bpoints, points, dim*(deg+1)*sizeof(float) );
      bufp += 2 + sizeof(short) + dim*(deg+1)*sizeof(float);
    }
    else if ( deg == 1 && N >0 ) {  /* a polyline */
      if ( bound[i].closing )
        bufp[0] = 1;
      else
        bufp[0] = 0;
      bufp[1] = (char)dim;
      *bn = (short)N;
      memcpy ( bpoints, points, dim*(N+1)*sizeof(float) );
      bufp += 2 + sizeof(short) + dim*(N+1)*sizeof(float);
    }
    else
      goto failure;
  }
  *bufp = 4;
  return buffer;

failure:
  pkv_SetScratchMemTop ( sp );
  return NULL;
} /*mbs_CompileTrimPatchBoundf*/

