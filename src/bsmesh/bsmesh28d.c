
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2014                            */
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

/* This procedure uses the procedures constructing the matrices of doubling */
/* and averaging corresponding to the direct doubling and averaging         */
/* procedures used by bsm_RefineBSMeshd. Therefore the correspondence       */
/* between the vertices of the input and output meshes is the same that the */
/* correspondence produced by the direct refinement procedure, which must   */
/* be preserved, should any modifications to this code be made.             */
boolean bsm_RefinementMatd ( int degree,
                             int inv, BSMvertex *imv, int *imvhei,
                             int inhe, BSMhalfedge *imhe,
                             int infac, BSMfacet *imfac, int *imfhei,
                             int *onv, BSMvertex **omv, int **omvhei,
                             int *onhe, BSMhalfedge **omhe,
                             int *onfac, BSMfacet **omfac, int **omfhei,
                             int *nrmat, index2 **rmi, double **rmc )
{
  void        *sp;
  int         nv1, nhe1, nfac1, *mvhei1, *mfhei1,
              nv2, nhe2, nfac2, *mvhei2, *mfhei2;
  BSMvertex   *mv1, *mv2;
  BSMhalfedge *mhe1, *mhe2;
  BSMfacet    *mfac1, *mfac2;
  index2      *mi1, *mi2, *mi3;
  double      *mc1, *mc2, *mc3;
  int         nnz1, nnz2, nnz3, *permut1, *permut2, *cols1, *cols2;
  int         i, nmult;

  sp = pkv_GetScratchMemTop ();

  mv1 = NULL;  mhe1 = NULL;  mfac1 = NULL;  mvhei1 = NULL;
  mfhei1 = NULL;  mi1 = NULL;  mc1 = NULL;
  mv2 = NULL;  mhe2 = NULL;  mfac2 = NULL;  mvhei2 = NULL;
  mfhei2 = NULL;  mi2 = NULL;  mc2 = NULL;
  mi3 = NULL;
  mc3 = NULL;

  if ( degree < 1 )
    goto failure;

  if ( !bsm_CheckMeshIntegrity ( inv, imv, imvhei, inhe, imhe,
                                 infac, imfac, imfhei, NULL, NULL ) )
    goto failure;

  bsm_DoublingNum ( inv, imv, imvhei, inhe, imhe, infac, imfac, imfhei,
                    &nv1, &nhe1, &nfac1 );
  nnz1 = bsm_DoublingMatSize ( inv, imv, imvhei, inhe, imhe, infac, imfac, imfhei );
  if ( nnz1 <= 0 )
    goto failure;
  PKV_MALLOC ( mv1, nv1*sizeof(BSMvertex) );
  PKV_MALLOC ( mhe1, nhe1*sizeof(BSMhalfedge) );
  PKV_MALLOC ( mfac1, nfac1*sizeof(BSMfacet) );
  PKV_MALLOC ( mvhei1, nhe1*sizeof(int) );
  PKV_MALLOC ( mfhei1, nhe1*sizeof(int) );
  PKV_MALLOC ( mi1, nnz1*sizeof(index2) );
  PKV_MALLOC ( mc1, nnz1*sizeof(double) );
  if ( !mv1 || !mhe1 || !mfac1 || !mvhei1 || !mfhei1 || !mi1 || !mc1 )
    goto failure;
  if ( !bsm_DoublingMatd ( inv, imv, imvhei, inhe, imhe, infac, imfac, imfhei,
                           &nv1, mv1, mvhei1, &nhe1, mhe1, &nfac1, mfac1, mfhei1,
                           &nnz1, mi1, mc1 ) )
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
    nnz2 = bsm_AveragingMatSize ( nv1, mv1, mvhei1, nhe1, mhe1,
                                  nfac1, mfac1, mfhei1 );
    if ( nnz2 <= 0 )
      goto failure;
    PKV_MALLOC ( mv2, nv2*sizeof(BSMvertex) );
    PKV_MALLOC ( mhe2, nhe2*sizeof(BSMhalfedge) );
    PKV_MALLOC ( mfac2, nfac2*sizeof(BSMfacet) );
    PKV_MALLOC ( mvhei2, nhe2*sizeof(int) );
    PKV_MALLOC ( mfhei2, nhe2*sizeof(int) );
    PKV_MALLOC ( mi2, nnz2*sizeof(index2) );
    PKV_MALLOC ( mc2, nnz2*sizeof(double) );
    if ( !mv2 || !mhe2 || !mfac2 || !mvhei2 || !mfhei2 || !mi2 || !mc2 )
      goto failure;
    if ( !bsm_AveragingMatd ( nv1, mv1, mvhei1, nhe1, mhe1, nfac1, mfac1, mfhei1,
                              &nv2, mv2, mvhei2, &nhe2, mhe2, &nfac2, mfac2, mfhei2,
                              &nnz2, mi2, mc2 ) )
      goto failure;
#ifdef DEBUG
{ boolean r;
  r = bsm_CheckMeshIntegrity ( nv2, mv2, mvhei2, nhe2, mhe2, nfac2, mfac2, mfhei2 );
  printf ( "averaging: %d\n", r );
  if ( !r )
    goto failure;
}
#endif
/* use the new, fast sparse matrix multiplication algorithm */
    permut1 = pkv_GetScratchMemi ( nnz1+nnz2+nv1+inv+2 );
    if ( !permut1 )
      goto failure;
    permut2 = &permut1[nnz1];
    cols1 = &permut2[nnz2];
    cols2 = &cols1[inv+1];
    if ( !pkn_SPMCountMMnnzC ( nv2, nv1, inv,
                               nnz2, mi2, (unsigned int*)permut2, cols2, false,
                               nnz1, mi1, (unsigned int*)permut1, cols1, false,
                               (unsigned int*)&nnz3, (unsigned int*)&nmult ) )
      goto failure;
    if ( nnz3 <= 0 )
      goto failure;

    PKV_MALLOC ( mi3, nnz3*sizeof(index2) );
    PKV_MALLOC ( mc3, nnz3*sizeof(double) );
    if ( !mi3 || !mc3 )
      goto failure;
    if ( !pkn_SPMmultMMCd ( nv2, nv1, inv,
                            nnz2, mi2, mc2, (unsigned int*)permut2, cols2, true,
                            nnz1, mi1, mc1, (unsigned int*)permut1, cols1, true,
                            mi3, mc3 ) )
      goto failure;
/* that's it */

    pkv_SetScratchMemTop ( sp );
    PKV_FREE ( mv1 );
    PKV_FREE ( mhe1 );
    PKV_FREE ( mfac1 );
    PKV_FREE ( mvhei1 );
    PKV_FREE ( mfhei1 );
    PKV_FREE ( mi1 );
    PKV_FREE ( mc1 );
    PKV_FREE ( mi2 );
    PKV_FREE ( mc2 );

    mv1 = mv2;        mv2 = NULL;
    mhe1 = mhe2;      mhe2 = NULL;
    mfac1 = mfac2;    mfac2 = NULL;
    mvhei1 = mvhei2;  mvhei2 = NULL;
    mfhei1 = mfhei2;  mfhei2 = NULL;
    nv1 = nv2;
    nhe1 = nhe2;
    nfac1 = nfac2;
    nnz1 = nnz3;
    mi1 = mi3;        mi3 = NULL;
    mc1 = mc3;        mc3 = NULL;
  }

  *onv = nv1;
  *onhe = nhe1;
  *onfac = nfac1;
  *omv = mv1;
  *omhe = mhe1;
  *omfac = mfac1;
  *omfhei = mfhei1;
  *omvhei = mvhei1;
  *nrmat = nnz1;
  *rmi = mi1;
  *rmc = mc1;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  *onv = *onhe = *onfac = 0;
  *omv = NULL;
  *omhe = NULL;
  *omfac = NULL;
  *omvhei = *omfhei = NULL;
  *nrmat = 0;
  *rmi = NULL;
  *rmc = NULL;
  if ( mv1 )     PKV_FREE ( mv1 );
  if ( mhe1 )    PKV_FREE ( mhe1 );
  if ( mfac1 )   PKV_FREE ( mfac1 );
  if ( mvhei1 )  PKV_FREE ( mvhei1 );
  if ( mfhei1 )  PKV_FREE ( mfhei1 );
  if ( mi1 )     PKV_FREE ( mi1 );
  if ( mc1 )     PKV_FREE ( mc1 );
  if ( mv2 )     PKV_FREE ( mv2 );
  if ( mhe2 )    PKV_FREE ( mhe2 );
  if ( mfac2 )   PKV_FREE ( mfac2 );
  if ( mvhei2 )  PKV_FREE ( mvhei2 );
  if ( mfhei2 )  PKV_FREE ( mfhei2 );
  if ( mi2 )     PKV_FREE ( mi2 );
  if ( mc2 )     PKV_FREE ( mc2 );
  if ( mi3 )     PKV_FREE ( mi3 );
  if ( mc3 )     PKV_FREE ( mc3 );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_RefinementMatd*/

