
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
boolean bsm_FindSpecialFSubnets ( int nv, BSMvertex *mv, int *mvhei,
                                  int nhe, BSMhalfedge *mhe,
                                  int nfac, BSMfacet *mfac, int *mfhei,
                                  int d, void *usrptr,
                                  void (*output)( int d, int k, int *vertnum,
                                                  int *mtab, void *usrptr ) )
{
#define MI(l,i,j) (((l)*(d+d+3)+(i))*(d+d+1)+(j))
  void *sp;
  int  i, j, k, l, m, p, fhe, he, f, v, vfhe, phe, ffhe;
  int  *mtab, *vertnum;

  if ( d < 1 )
    return true;
  sp = pkv_GetScratchMemTop ();
  for ( i = 0; i < nfac; i++ ) {
    k = mfac[i].degree;
    if ( k == 4 )
      goto nextfac;  /* process special facets only */
    mtab = pkv_GetScratchMemi ( k*(d+d+1)*(d+d+3) );
    vertnum = pkv_GetScratchMemi ( k*(d+1)*(d+1) );
    if ( !mtab || !vertnum )
      goto failure;
memset ( mtab, -1, k*(d+d+1)*(d+d+3)*sizeof(int) );
    fhe = mfac[i].firsthalfedge;
          /* find k regular subnets */
    for ( l = 0; l < k; l++ ) {
      mtab[MI(l,1,0)] = he = mhe[mfhei[fhe+l]].otherhalf;
      if ( he < 0 )
        goto nextfac;
      mtab[MI(l,2,0)] = mhe[he].v0;
      mtab[MI(l,0,0)] = mhe[he].v1;
      mtab[MI(l,1,1)] = f = mhe[he].facetnum;
      if ( mfac[f].degree != 4 )
        goto nextfac;
      ffhe = mfac[f].firsthalfedge;
      for ( p = 0; p < 4; p++ )
        if ( mfhei[ffhe+p] == he )
          break;
      if ( p >= 4 )
        goto failure;
      mtab[MI(l,0,1)] = mfhei[ffhe+(p+1) % 4];
      mtab[MI(l,1,2)] = he = mfhei[ffhe+(p+2) % 4];
      mtab[MI(l,0,2)] = mhe[he].v0;
      mtab[MI(l,2,1)] = he = mfhei[ffhe+(p+3) % 4];
      mtab[MI(l,2,2)] = mhe[he].v0;
            /* go north */
      for ( j = 3; j <= d+d; j += 2 ) {
        v = mtab[MI(l,0,j-1)];
        if ( mv[v].degree != 4 )
          goto nextfac;
        vfhe = mv[v].firsthalfedge;
        if ( mhe[mvhei[vfhe+3]].otherhalf < 0 )
          goto nextfac;
        phe = mhe[mtab[MI(l,0,j-2)]].otherhalf;
        for ( p = 0; p < 4; p++ )
          if ( mvhei[vfhe+p] == phe )
            break;
        if ( p >= 4 )
          goto failure;
        mtab[MI(l,0,j)] = he = mvhei[vfhe+(p+2) % 4];
        mtab[MI(l,1,j)] = f = mhe[he].facetnum;
        if ( mfac[f].degree != 4 )
          goto nextfac;
        mtab[MI(l,0,j+1)] = mhe[he].v1;
        ffhe = mfac[f].firsthalfedge;
        he = mhe[mtab[MI(l,1,j-1)]].otherhalf;
        for ( p = 0; p < 4; p++ )
          if ( mfhei[ffhe+p] == he )
            break;
        if ( p >= 4 )
          goto failure;
        mtab[MI(l,2,j)] = he = mfhei[ffhe+(p+3) % 4];
        mtab[MI(l,2,j+1)] = mhe[he].v0;
        mtab[MI(l,1,j+1)] = mfhei[ffhe+(p+2) % 4];
      }
            /* go east */
      for ( m = 3; m <= d+d+2; m += 2 ) {
        v = mtab[MI(l,m-1,0)];
        if ( mv[v].degree != 4 )
          goto nextfac;
        vfhe = mv[v].firsthalfedge;
        if ( mhe[mvhei[vfhe+3]].otherhalf < 0 )
          goto nextfac;
        phe = mhe[mtab[MI(l,m-1,1)]].otherhalf;
        for ( p = 0; p < 4; p++ )
          if ( mvhei[vfhe+p] == phe )
            break;
        if ( p >= 4 )
          goto failure;
        mtab[MI(l,m,0)] = he = mhe[mvhei[vfhe+(p+3) % 4]].otherhalf;
        mtab[MI(l,m+1,0)] = mhe[he].v0;
        for ( j = 1; j <= d+d; j += 2 ) {
          phe = mtab[MI(l,m,j-1)];
          v = mhe[phe].v1;
          if ( mv[v].degree != 4 )
            goto nextfac;
          mtab[MI(l,m,j)] = f = mhe[phe].facetnum;
          if ( mfac[f].degree != 4 )
            goto nextfac;
          ffhe = mfac[f].firsthalfedge;
          for ( p = 0; p < 4; p++ )
            if ( mfhei[ffhe+p] == phe )
              break;
          if ( p >= 4 )
            goto failure;
          mtab[MI(l,m+1,j)] = he = mfhei[ffhe+(p+3) % 4];
          mtab[MI(l,m+1,j+1)] = mhe[he].v0;
          mtab[MI(l,m,j+1)] = mhe[mfhei[ffhe+(p+2) % 4]].otherhalf;
        }
      }
    }
            /* verify, whether the regular subnets fit */
    for ( l = 0, m = k-1;  l < k;  m = l++ ) {
      for ( j = 1; j < d+d; j += 2 )
        if ( mtab[MI(m,j+2,0)] != mhe[mtab[MI(l,0,j)]].otherhalf )
          goto nextfac;
    }
            /* extract vertex numbers and output */
    for ( l = m = 0;  l < k;  l++ )
      for ( j = 0; j <= d+d; j += 2 )
        for ( p = 2;  p <= d+d+2;  p += 2, m ++ )
          vertnum[m] = mtab[MI(l,p,j)];
    output ( d, k, vertnum, mtab, usrptr );
nextfac:
    pkv_SetScratchMemTop ( sp );
  }
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef MI
} /*bsm_FindSpecialFSubnets*/

