
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
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

/*#define DEBUG*/

/* This procedure uses the procedures of direct doubling and averaging,     */
/* which produce the same correspondence between the vertices of the input  */
/* and output meshes that the procedures constructing the matrices of these */
/* transformations, used by bsm_RefinementMatd. This correspondence must be */
/* preserved, should any modifications to this code be made.             */
boolean bsm_RefineBSMeshd ( int spdimen, int degree,
                    int inv, BSMvertex *imv, int *imvhei, double *iptc,
                    int inhe, BSMhalfedge *imhe,
                    int infac, BSMfacet *imfac, int *imfhei,
                    int *onv, BSMvertex **omv, int **omvhei, double **optc,
                    int *onhe, BSMhalfedge **omhe,
                    int *onfac, BSMfacet **omfac, int **omfhei )
{
  void        *sp;
  int         nv1, nhe1, nfac1, *mvhei1, *mfhei1,
              nv2, nhe2, nfac2, *mvhei2, *mfhei2;
  BSMvertex   *mv1, *mv2;
  BSMhalfedge *mhe1, *mhe2;
  BSMfacet    *mfac1, *mfac2;
  double      *ptc1, *ptc2;
  int         i;

  sp = pkv_GetScratchMemTop ();

  mv1 = NULL;  mhe1 = NULL;  mfac1 = NULL;  mvhei1 = NULL;
  mfhei1 = NULL;  ptc1 = NULL;
  mv2 = NULL;  mhe2 = NULL;  mfac2 = NULL;  mvhei2 = NULL;
  mfhei2 = NULL;  ptc2 = NULL;

  if ( spdimen < 0 || degree < 1 )
    goto failure;

  if ( !bsm_CheckMeshIntegrity ( inv, imv, imvhei, inhe, imhe,
                                 infac, imfac, imfhei ) )
    goto failure;

  bsm_DoublingNum ( inv, imv, imvhei, inhe, imhe, infac, imfac, imfhei,
                    &nv1, &nhe1, &nfac1 );
  PKV_MALLOC ( mv1, nv1*sizeof(BSMvertex) );
  PKV_MALLOC ( mhe1, nhe1*sizeof(BSMhalfedge) );
  PKV_MALLOC ( mfac1, nfac1*sizeof(BSMfacet) );
  PKV_MALLOC ( mvhei1, nhe1*sizeof(int) );
  PKV_MALLOC ( mfhei1, nhe1*sizeof(int) );
  PKV_MALLOC ( ptc1, spdimen*nv1*sizeof(double) );
  if ( !mv1 || !mhe1 || !mfac1 || !mvhei1 || !mfhei1 || !ptc1 )
    goto failure;
  if ( !bsm_Doublingd ( spdimen, inv, imv, imvhei, iptc,
                        inhe, imhe, infac, imfac, imfhei,
                        &nv1, mv1, mvhei1, ptc1, &nhe1, mhe1,
                        &nfac1, mfac1, mfhei1 ) )
    goto failure;
#ifdef DEBUG
{ boolean r;
  r = bsm_CheckMeshIntegrity ( nv1, mv1, mvhei1, nhe1, mhe1, nfac1, mfac1, mfhei1 );
  printf ( "doubling: %d\n", r );
  if ( !r )
    goto failure;
}
#endif

  for ( i = 0; i < degree; i++ ) {
    bsm_AveragingNum ( nv1, mv1, mvhei1, nhe1, mhe1, nfac1, mfac1, mfhei1,
                       &nv2, &nhe2, &nfac2 );
    if ( nv2 <= 0 || nhe2 <= 0 || nfac2 <= 0 )
      goto failure;
    PKV_MALLOC ( mv2, nv2*sizeof(BSMvertex) );
    PKV_MALLOC ( mhe2, nhe2*sizeof(BSMhalfedge) );
    PKV_MALLOC ( mfac2, nfac2*sizeof(BSMfacet) );
    PKV_MALLOC ( mvhei2, nhe2*sizeof(int) );
    PKV_MALLOC ( mfhei2, nhe2*sizeof(int) );
    PKV_MALLOC ( ptc2, spdimen*nv2*sizeof(double) );
    if ( !mv2 || !mhe2 || !mfac2 || !mvhei2 || !mfhei2 || !ptc2 )
      goto failure;
    if ( !bsm_Averagingd ( spdimen, nv1, mv1, mvhei1, ptc1,
                           nhe1, mhe1, nfac1, mfac1, mfhei1,
                           &nv2, mv2, mvhei2, ptc2, &nhe2, mhe2,
                           &nfac2, mfac2, mfhei2 ) )
      goto failure;
#ifdef DEBUG
{ boolean r;
  r = bsm_CheckMeshIntegrity ( nv2, mv2, mvhei2, nhe2, mhe2, nfac2, mfac2, mfhei2 );
  printf ( "averaging: %d\n", r );
  if ( !r )
    goto failure;
}
#endif
    PKV_FREE ( mv1 );
    PKV_FREE ( mhe1 );
    PKV_FREE ( mfac1 );
    PKV_FREE ( mvhei1 );
    PKV_FREE ( mfhei1 );
    PKV_FREE ( ptc1 );
    mv1 = mv2;        mv2 = NULL;
    mhe1 = mhe2;      mhe2 = NULL;
    mfac1 = mfac2;    mfac2 = NULL;
    mvhei1 = mvhei2;  mvhei2 = NULL;
    mfhei1 = mfhei2;  mfhei2 = NULL;
    ptc1 = ptc2;      ptc2 = NULL;
    nv1 = nv2;
    nhe1 = nhe2;
    nfac1 = nfac2;
  }

  *onv = nv1;
  *onhe = nhe1;
  *onfac = nfac1;
  *omv = mv1;
  *omhe = mhe1;
  *omfac = mfac1;
  *omfhei = mfhei1;
  *omvhei = mvhei1;
  *optc = ptc1;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  *onv = *onhe = *onfac = 0;
  *omv = NULL;
  *omhe = NULL;
  *omfac = NULL;
  *omvhei = *omfhei = NULL;
  *optc = NULL;
  if ( mv1 )    PKV_FREE ( mv1 );
  if ( mhe1 )   PKV_FREE ( mhe1 );
  if ( mfac1 )  PKV_FREE ( mfac1 );
  if ( mvhei1 ) PKV_FREE ( mvhei1 );
  if ( mfhei1 ) PKV_FREE ( mfhei1 );
  if ( ptc1 )   PKV_FREE ( ptc1 );
  if ( mv2 )    PKV_FREE ( mv2 );
  if ( mhe2 )   PKV_FREE ( mhe2 );
  if ( mfac2 )  PKV_FREE ( mfac2 );
  if ( mvhei2 ) PKV_FREE ( mvhei2 );
  if ( mfhei2 ) PKV_FREE ( mfhei2 );
  if ( ptc2 )   PKV_FREE ( ptc2 );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_RefineBSMeshd*/

