
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>   /* for testing */

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"

#undef CONST_
#define CONST_

#include "multibs.h"

#define WGT 0.5*SQRT2

/* ///////////////////////////////////////////////////////////////////////// */
/* It is assumed that the array outknots contains elements */
/* (knots) on positions from -outdeg to Nout+outdeg.       */
static boolean _mbs_DegRedFindMatrixf ( int indegree,
                     int Naux, float *auxknots, int deltadeg,
                     int outdeg, int Nout, float *outknots,
                     int *nrows, int *ncols,
                     bandm_profile **prof, float **mtrx )
{
  void          *sp;
  int           i, j, _nrows, _ncols, coll, fkn, kmat;
  int           ii, jj, kk, ll, mm, ra, rb;
  float         ta, tb;
  bandm_profile *_prof;
  float         *_mtrx, *a, *b, *auxknots2;

  *nrows = _nrows = Naux-indegree;
  *ncols = _ncols = Nout-outdeg;
  coll = deltadeg*(outdeg+1)+1;
  *prof = _prof = pkv_GetScratchMem ( (_ncols+1)*sizeof(bandm_profile) );
  *mtrx = _mtrx = pkv_GetScratchMemf ( _ncols*coll );
  sp = pkv_GetScratchMemTop ();
  a = pkv_GetScratchMemf ( 2*outdeg+1 );
  b = pkv_GetScratchMemf ( Naux );
  auxknots2 = pkv_GetScratchMemf ( 3*outdeg*(deltadeg+1)+2 );
  if ( !_prof || !_mtrx || !a || !b || !auxknots2 ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    return false;
  }

  memset ( a, 0, (2*outdeg+1)*sizeof(float) );
  a[outdeg] = 1.0;
  for ( i = fkn = kmat = 0; i < _ncols; i++ ) {
        /* compute the matrix column coefficients */
    mbs_BSDegElevC1f ( outdeg, 3*outdeg+1, &outknots[i-outdeg], a, deltadeg,
                       &ii, &mm, auxknots2, b, true );
        /* identify the nonzero ones and their position in the column */
    for ( j = i; auxknots[0] > outknots[j]; j++ ) ;
    ta = outknots[j];  tb = outknots[i+outdeg+1];
    for ( ra = 0;  outknots[j+ra+1] == ta;  ra++ ) ;
    for ( rb = 0;  outknots[i+outdeg-rb] == tb;  rb++ ) ;
    for ( ii = deltadeg; auxknots2[ii+1] <= ta;  ii++ ) ;
    for ( jj = mm-deltadeg; auxknots2[jj-1] >= tb;  jj-- ) ;
    kk = ii-ra-deltadeg ;  ll = jj+rb-outdeg-1;
    for ( mm = fkn; auxknots[mm+1] <= ta; mm++ ) ;
    fkn = mm-ii+kk;
    if ( fkn < 0 ) { kk -= fkn;  fkn = 0; }
    if ( fkn+ll-kk+1 >= _nrows ) ll = _nrows-fkn+kk-1;
        /* store them in the band matrix representation */
    _prof[i].firstnz = fkn;
    _prof[i].ind = kmat;
    memcpy ( &_mtrx[kmat], &b[kk], (ll-kk+1)*sizeof(float) );
    kmat += ll-kk+1;
  }
  _prof[_ncols].firstnz = 0;
  _prof[_ncols].ind = kmat;
  pkv_SetScratchMemTop ( sp );
  return true;
} /*_mbs_DegRedFindMatrixf*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean _mbs_FindDegRedKnotSequencesf ( int indeg, int inlastknot,
                 const float *inknots, int deltadeg,
                 int *auxlastknot, float **auxknots,
                 int *outdegree, int *outlastknot, float **outknots )
{
  int   outdeg, Nout, Naux;
  int   i, j, k, mult;
  float t, *auxkn, *outkn;

  *outdegree = outdeg = indeg-deltadeg;

        /* compute the lengths of the auxiliary and final knot sequences */
  Naux = Nout = 0;
  for ( i = 0; i <= inlastknot; i += mult ) {
    t = inknots[i];  mult = 1;
    while ( i+mult <= inlastknot && t == inknots[i+mult] )
      mult ++;
    j = max ( mult, deltadeg+1 );
    Naux += j;
    j = min ( j-deltadeg, outdeg+1 );
    Nout += j;
  }

        /* allocate the array for the knot sequences,  */
        /* leaving space for knot sequence corrections */
  *auxknots = auxkn = pkv_GetScratchMemf ( Naux );
  outkn = pkv_GetScratchMemf ( Nout+4*outdeg );
  if ( !auxkn || !outkn )
    return false;

        /* generate the knot sequences */
  for ( i = 0; i < outdeg; i++ )
    outkn[i] = inknots[0];
  *outknots = outkn = &outkn[outdeg];
  Naux = Nout = 0;
  for ( i = 0; i <= inlastknot; i += mult ) {
    t = inknots[i];  mult = 1;
    while ( i+mult <= inlastknot && t == inknots[i+mult] )
      mult ++;
    j = max ( mult, deltadeg+1 );
    for ( k = 0; k < j; k++ )
      auxkn[Naux++] = t;
    j = min ( j-deltadeg, outdeg+1 );
    for ( k = 0; k < j; k++ )
      outkn[Nout++] = t;
  }
        /* make corrections if necessary */
  for ( t = inknots[indeg], i = 0;  outkn[outdeg-i] > t;  i++ ) ;
  if ( i ) {
    memmove ( &outkn[i], outkn, Nout*sizeof(float) );  Nout += i;
    for ( t = inknots[0]; i; outkn[--i] = t ) ;
  }
  for ( t = inknots[indeg], i = 0; outkn[outdeg+i+1] == t; i++ );
  if ( i ) {
    Nout -= i;
    memmove ( outkn, &outkn[i], Nout*sizeof(float) );
  }
  for ( t = inknots[inlastknot-indeg], i = 0;  outkn[Nout-outdeg+i] < t;  i++ ) ;
  if ( i )
    for ( t = inknots[inlastknot]; i; outkn[++Nout] = t, i-- ) ;
  for ( t = inknots[inlastknot-indeg];  outkn[Nout-outdeg-2] >= t;  Nout-- ) ;

  for ( i = 0; i < outdeg; i++ )
    outkn[Nout+i] = inknots[inlastknot];

  *auxlastknot = Naux-1;
  *outlastknot = Nout-1;
  return true;
} /*_mbs_FindDegRedKnotSequencesf*/

boolean mbs_multiBSDegRedf ( int ncurves, int spdimen,
                             int indegree, int inlastknot, const float *inknots,
                             int inpitch, CONST_ float *inctlpoints,
                             int deltadeg,
                             int *outdegree, int *outlastknot, float *outknots,
                             int outpitch, float *outctlpoints )
{
  void          *sp;
  int           i, j, Naux, Nout, outdeg, auxpitch, nrows, ncols;
  float         *auxknots, *outkn, *mtrx, *a;
  bandm_profile *prof;

  sp = pkv_GetScratchMemTop ();
  if ( deltadeg <= 0 || deltadeg > indegree )
    goto failure;

  if ( !_mbs_FindDegRedKnotSequencesf ( indegree, inlastknot, inknots,
             deltadeg, &Naux, &auxknots, &outdeg, &Nout, &outkn ) )
    goto failure;

  auxpitch = (Naux-indegree)*spdimen;
  a = pkv_GetScratchMemf ( ncurves*auxpitch );
  if ( !a )
    goto failure;
  mbs_multiOsloInsertKnotsf ( ncurves, spdimen, indegree, inlastknot, inknots,
                          inpitch, inctlpoints, Naux, auxknots, auxpitch, a );
  mbs_multiRemoveSuperfluousKnotsf ( ncurves, spdimen, indegree,
                                     &Naux, auxknots, auxpitch, auxpitch, a );

        /* skip the degenerate knot intervals at the domain ends */
  for ( i = 1; outkn[outdeg+i] == outkn[outdeg]; i++ ) ;
  for ( j = 1; auxknots[indegree+j] <= outkn[outdeg]; j++ ) ;
  if ( (j -= i) > 0 )
    { auxknots += j;  Naux -= j; a = &a[j*spdimen]; }
  for ( i = 1; outkn[Nout-outdeg-i] == outkn[Nout-outdeg]; i++ ) ;
  for ( j = 1; auxknots[Naux-indegree-j] >= outkn[Nout-outdeg]; j++ ) ;
  if ( (j -= i) > 0 )
    Naux -= j;

  _mbs_DegRedFindMatrixf ( indegree, Naux, auxknots, deltadeg, outdeg,
                           Nout, outkn, &nrows, &ncols, &prof, &mtrx );

      /* solve the least-squares problem to obtain the resulting control points */
  pkn_multiBandmSolveRLSQf ( nrows, ncols, prof, mtrx,
                             ncurves, spdimen, auxpitch, a,
                             outpitch, outctlpoints );

  memcpy ( outknots, outkn, (Nout+1)*sizeof(float) );
  *outdegree = outdeg;
  *outlastknot = Nout;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiBSDegRedf*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean _mbs_FindClosedDegRedKnotSequencesf ( int indeg, int inlastknot,
                 const float *inknots, int deltadeg,
                 int *auxlastknot, float **auxknots,
                 int *outdegree, int *outlastknot, float **outknots,
                 int *kaux, int *kout )
{
  int   outdeg, Nout, Naux;
  int   i, j, k, mult, Kin, Kaux, Kout;
  float t, T, *auxkn, *outkn;

  *outdegree = outdeg = indeg-deltadeg;
  Kin = inlastknot-2*indeg;
  T = inknots[indeg+Kin]-inknots[indeg];

        /* compute the lengths of the auxiliary and final knot sequences */
  Naux = Nout = 0;
  for ( i = 0; i <= inlastknot; i += mult ) {
    t = inknots[i];  mult = 1;
    while ( i+mult <= inlastknot && t == inknots[i+mult] )
      mult ++;
    j = max ( mult, deltadeg+1 );
    Naux += j;
    j = min ( j-deltadeg, outdeg+1 );
    Nout += j;
  }

        /* allocate the array for the knot sequences,  */
        /* leaving space for knot sequence corrections */
  *auxknots = auxkn = pkv_GetScratchMemf ( Naux );
  outkn = pkv_GetScratchMemf ( Nout+4*outdeg );
  if ( !auxkn || !outkn )
    return false;

        /* generate the knot sequences */
  *outknots = outkn = &outkn[outdeg];
  Naux = Nout = 0;
  for ( i = 0; i <= inlastknot; i += mult ) {
    t = inknots[i];  mult = 1;
    while ( i+mult <= inlastknot && t == inknots[i+mult] )
      mult ++;
    j = max ( mult, deltadeg+1 );
    for ( k = 0; k < j; k++ )
      auxkn[Naux++] = t;
    j = min ( j-deltadeg, outdeg+1 );
    for ( k = 0; k < j; k++ )
      outkn[Nout++] = t;
  }
  for ( i = 0; auxkn[i+1] <= inknots[indeg]; i++ ) ;
  for ( j = Naux; auxkn[j-1] > inknots[inlastknot-indeg]; j-- ) ;
  Kaux = j-i-1;
  for ( i = 0; outkn[i+1] <= inknots[indeg]; i++ ) ;
  for ( j = Nout; outkn[j-1] > inknots[inlastknot-indeg]; j-- ) ;
  Kout = j-i-1;

        /* make corrections if necessary */
  for ( t = inknots[indeg], i = 0;  outkn[outdeg-i] > t;  i++ ) ;
  if ( i )
    { memmove ( &outkn[i], outkn, Nout*sizeof(float) );  Nout += i; }
  else {
    for ( t = inknots[indeg], i = 0;  outkn[outdeg+i+1] == t;  i++ ) ;
    if ( i ) {
      Nout -= i;
      memmove ( outkn, &outkn[i], Nout*sizeof(float) );
    }
    i = 0;
  }
  for ( i--; i >= -outdeg; i-- ) {
    t = outkn[i+Kout] - T;
    outkn[i] = min ( t, outkn[i+1] );  /* for rounding errors */
  }
  for ( t = inknots[inlastknot-indeg], i = 0;  outkn[Nout-outdeg+i] < t;  i++ ) ;
  if ( i )
    for ( i--; i >= 0; i-- ) {
      t = outkn[Nout-Kout] + T;
      outkn[Nout] = max ( t, outkn[Nout-1] );  /* for rounding errors */
      Nout ++;
    }
  while ( Nout <= Kout+3*outdeg ) {
    t = outkn[Nout-Kout] + T;
    outkn[Nout] = max ( t, outkn[Nout-1] );  /* for rounding errors */
    Nout++;
  }

  *auxlastknot = Naux-1;
  *outlastknot = Nout-1;
  *kaux = Kaux;
  *kout = Kout;
  return true;
} /*_mbs_FindClosedDegRedKnotSequencesf*/

boolean mbs_multiBSDegRedClosedf ( int ncurves, int spdimen,
                             int indegree, int inlastknot, const float *inknots,
                             int inpitch, CONST_ float *inctlpoints,
                             int deltadeg,
                             int *outdegree, int *outlastknot, float *outknots,
                             int outpitch, float *outctlpoints )
{
  void          *sp;
  int           Naux, Nout, outdeg, auxpitch, nrows, ncols, Kaux, Kout;
  int           i, j, k, l, m, p, q;
  float         *auxknots, *outkn, *mtrx, *a, *c;
  bandm_profile *prof;

  sp = pkv_GetScratchMemTop ();
  if ( deltadeg <= 0 || deltadeg > indegree )
    goto failure;

  if ( !_mbs_FindClosedDegRedKnotSequencesf ( indegree, inlastknot, inknots,
          deltadeg, &Naux, &auxknots, &outdeg, &Nout, &outkn, &Kaux, &Kout ) )
    goto failure;

  auxpitch = (Naux-indegree)*spdimen;
  a = pkv_GetScratchMemf ( ncurves*auxpitch );
  if ( !a )
    goto failure;
  mbs_multiOsloInsertKnotsf ( ncurves, spdimen, indegree, inlastknot, inknots,
                          inpitch, inctlpoints, Naux, auxknots, auxpitch, a );
  mbs_multiRemoveSuperfluousKnotsf ( ncurves, spdimen, indegree,
                                     &Naux, auxknots, auxpitch, auxpitch, a );
        /* skip the degenerate knot intervals at the domain ends */
  for ( i = 1; outkn[outdeg+i] == outkn[outdeg]; i++ ) ;
  for ( j = 1; auxknots[indegree+j] <= outkn[outdeg]; j++ ) ;
  if ( (j -= i) > 0 )
    { auxknots += j;  a = &a[j*spdimen]; }
  Nout = Kout+2*outdeg;
  Naux = Kaux+2*indegree;
  for ( j = 1; auxknots[Naux-indegree-j] >= inknots[inlastknot-indegree]; j++ ) ;
  if ( (p = j - i) > 0 )
    Naux -= p;
  for ( q = 0; outkn[Nout-outdeg-q-1] == outkn[Nout-outdeg]; q++ ) ;

  _mbs_DegRedFindMatrixf ( indegree, Naux, auxknots, deltadeg, outdeg,
                           Nout-q, outkn, &nrows, &ncols, &prof, &mtrx );

        /* apply weights */
  pkn_MultMatrixNumf ( ncurves, (indegree-p)*spdimen, auxpitch,
                       a, WGT, auxpitch, a );
  pkn_MultMatrixNumf ( ncurves, (indegree-p)*spdimen, auxpitch,
                       &a[Kaux*spdimen], WGT, auxpitch, &a[Kaux*spdimen] );
  for ( i = 0; i < ncols; i++ ) {
    k = prof[i].firstnz;
    l = prof[i].firstnz + prof[i+1].ind - prof[i].ind;
    for ( j = k, m = min(l,indegree-p);  j < m;  j++ )
      mtrx[prof[i].ind+j-k] *= (float)(WGT);
    for ( j = min (Naux-indegree-1,l-1), m = max(k,Kaux);  j >= m;  j-- )
      mtrx[prof[i].ind+j-k] *= (float)(WGT);
  }

      /* setup the constraints matrix */
  if ( !(c = pkv_GetScratchMemf ( (outdeg-q)*(Nout-outdeg-q) )) ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  memset ( c, 0, (outdeg-q)*(Nout-outdeg-q)*sizeof(float) );
  for ( i = 0; i < outdeg-q; i++ ) {
    c[i*(Nout-outdeg-q+1)] = 1.0;
    c[i*(Nout-outdeg-q+1)+Kout] = -1.0;
  }

      /* solve the linear least-squares problem with constraints */
  pkn_multiBandmSolveCRLSQf ( nrows, ncols, prof, mtrx, outdeg-q,
                              Nout-outdeg-q, c, ncurves, spdimen, auxpitch, a,
                              0, NULL, outpitch, outctlpoints );
  if ( q > 0 )
    pkv_Movef ( ncurves, q*spdimen, outpitch, Kout*spdimen,
                &outctlpoints[(outdeg-q)*spdimen] );

  memcpy ( outknots, outkn, (Nout+1)*sizeof(float) );
  *outdegree = outdeg;
  *outlastknot = Nout;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiBSDegRedClosedf*/

