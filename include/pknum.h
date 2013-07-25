
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libpknum library of C procedures - numerical methods */

#ifndef PKNUM_H
#define PKNUM_H

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ///////////////////////////////////////////////////////////////////////// */

typedef struct bandm_profile {
    int firstnz;  /* index of the first row with a nonzero  */
                  /* coefficient in the column */
    int ind;      /* index of this coefficient in the array */
  } bandm_profile;

typedef struct {  /* structure for representing sparse matrices */
    int i, j;
  } index2;

typedef struct {  /* structure for representing sparse submatrices */
    int i, j, k;
  } index3;

typedef struct { float x, y; } complexf;
typedef struct { double x, y; } complexd;


#define pkn_LowerTrMatIndex(i,j) \
  ( (i)*((i)+1)/2+(j) )
#define pkn_SymMatIndex(i,j) \
  ( (i) >= (j) ? (i)*((i)+1)/2+(j) : (j)*((j)+1)/2+(i) )
#define pkn_LHessenbergMatIndex(i,j) \
  (pkn_LowerTrMatIndex(((i)+1),(j))-1)
#define pkn_UHessenbergMatIndex(i,j) \
  (pkn_LowerTrMatIndex(((j)+1),(i))-1)

#include "pknumf.h"
#include "pknumd.h"

/* ///////////////////////////////////////////////////////////////////////// */
void pkn_BandmFindQRMSizes ( int ncols, const bandm_profile *aprof,
                             int *qsize, int *rsize );

void pkn_PrintProfile ( int ncols, const bandm_profile *prof );

/* ///////////////////////////////////////////////////////////////////////// */
int pkn_Block1ArraySize ( int k, int r, int s );
int pkn_Block1FindBlockPos ( int k, int r, int s, int i, int j );
int pkn_Block1FindElemPos ( int k, int r, int s, int i, int j );

/* ///////////////////////////////////////////////////////////////////////// */
int pkn_Block2ArraySize ( int k, int r, int s, int t );
int pkn_Block2FindBlockPos ( int k, int r, int s, int t, int i, int j );
int pkn_Block2FindElemPos ( int k, int r, int s, int t, int i, int j );

/* ///////////////////////////////////////////////////////////////////////// */
int pkn_Block3ArraySize ( int k, int r, int s );
int pkn_Block3FindBlockPos ( int k, int r, int s, int i, int j );
int pkn_Block3FindElemPos ( int k, int r, int s, int i, int j ); 

/* ///////////////////////////////////////////////////////////////////////// */
int pkn_NRBArraySize ( int n, const int *prof );

/* ///////////////////////////////////////////////////////////////////////// */
int pkn_TMBSize ( int n );
boolean pkn_TMBElem ( byte *bittm, int i, int j );
void pkn_TMBElemSet ( byte *bittm, int i, int j );
void pkn_TMBElemClear ( byte *bittm, int i, int j );
boolean pkn_TMBTestAndSet ( byte *bittm, int i, int j );
boolean pkn_TMBTestAndClear ( byte *bittm, int i, int j );

/* ///////////////////////////////////////////////////////////////////////// */
/* auxiliary procedures */
void pkn_SPMindex2to3 ( unsigned int nnz, index2 *ai, index3 *sai );
void pkn_SPMindex3to2 ( unsigned int nnz, index3 *sai, index2 *ai );

boolean pkn_SPMSortByRows ( unsigned int nrows, unsigned int ncols,
                            unsigned int nnz, index2 *ai, unsigned int *permut );
boolean pkn_SPMSortByCols ( unsigned int nrows, unsigned int ncols,
                            unsigned int nnz, index2 *ai, unsigned int *permut );
boolean pkn_SPMFindRows ( unsigned int nrows, unsigned int ncols, unsigned int nnz,
                          index2 *ai, unsigned int *permut, boolean ro,
                          int *rows );
boolean pkn_SPMFindCols ( unsigned int nrows, unsigned int ncols, unsigned int nnz,
                          index2 *ai, unsigned int *permut, boolean co,
                          int *cols );

/* fast multiplication algorithm - uses considerable additional memory */
boolean pkn_SPMCountMMnnzR ( int nra, int nca, int ncb,
                             unsigned int nnza, index2 *ai,
                             unsigned int *apermut, int *arows, boolean ra,
                             unsigned int nnzb, index2 *bi,
                             unsigned int *bpermut, int *brows, boolean rb,
                             unsigned int *nnzab, unsigned int *nmultab );
boolean pkn_SPMFindMMnnzR ( int nra, int nca, int ncb,
                            int unsigned nnza, index2 *ai,
                            unsigned int *apermut, int *arows,
                            unsigned int nnzb, index2 *bi,
                            unsigned int *bpermut, int *brows,
                            index2 *abi, int *abpos, index2 *aikbkj );
boolean pkn_SPMCountMMnnzC ( int nra, int nca, int ncb,
                             unsigned int nnza, index2 *ai,
                             unsigned int *apermut, int *acols, boolean ca,
                             unsigned int nnzb, index2 *bi,
                             unsigned int *bpermut, int *bcols, boolean cb,
                             unsigned int *nnzab, unsigned int *nmultab );
boolean pkn_SPMFindMMnnzC ( int nra, int nca, int ncb,
                            unsigned int nnza, index2 *ai,
                            unsigned int *apermut, int *acols,
                            unsigned int nnzb, index2 *bi,
                            unsigned int *bpermut, int *bcols,
                            index2 *abi, int *abpos, index2 *aikbkj );

boolean pkn_SPMCountMMTnnzR ( int nra, int nca, int nrb,
                              unsigned int nnza, index2 *ai,
                              unsigned int *apermut, int *arows, boolean ra,
                              unsigned int nnzb, index2 *bi,
                              unsigned int *bpermut, int *bcols, boolean cb,
                              unsigned int *nnzab, unsigned int *nmultab );
boolean pkn_SPMFindMMTnnzR ( int nra, int nca, int nrb,
                             unsigned int nnza, index2 *ai,
                             unsigned int *apermut, int *arows,
                             unsigned int nnzb, index2 *bi,
                             unsigned int *bpermut, int *bcols,
                             index2 *abi, int *abpos, index2 *aikbkj );
boolean pkn_SPMCountMMTnnzC ( int nra, int nca, int nrb,
                              unsigned int nnza, index2 *ai,
                              unsigned int *apermut, int *acols, boolean ca,
                              unsigned int nnzb, index2 *bi,
                              unsigned int *bpermut, int *brows, boolean rb,
                              unsigned int *nnzab, unsigned int *nmultab );
boolean pkn_SPMFindMMTnnzC ( int nra, int nca, int nrb,
                             unsigned int nnza, index2 *ai,
                             unsigned int *apermut, int *acols,
                             unsigned int nnzb, index2 *bi,
                             unsigned int *bpermut, int *brows,
                             index2 *abi, int *abpos, index2 *aikbkj );

boolean pkn_SPMCountMTMnnzR ( int nra, int nca, int ncb,
                              unsigned int nnza, index2 *ai,
                              unsigned int *apermut, int *acols, boolean ca,
                              unsigned int nnzb, index2 *bi,
                              unsigned int *bpermut, int *brows, boolean rb,
                              unsigned int *nnzab, unsigned int *nmultab );
boolean pkn_SPMFindMTMnnzR ( int nra, int nca, int ncb,
                             unsigned int nnza, index2 *ai,
                             unsigned int *apermut, int *acols,
                             unsigned int nnzb, index2 *bi,
                             unsigned int *bpermut, int *brows,
                             index2 *abi, int *abpos, index2 *aikbkj );
boolean pkn_SPMCountMTMnnzC ( int nra, int nca, int ncb,
                              unsigned int nnza, index2 *ai,
                              unsigned int *apermut, int *arows, boolean ra,
                              unsigned int nnzb, index2 *bi,
                              unsigned int *bpermut, int *bcols, boolean cb,
                              unsigned int *nnzab, unsigned int *nmultab );
boolean pkn_SPMFindMTMnnzC ( int nra, int nca, int ncb,
                             unsigned int nnza, index2 *ai,
                             unsigned int *apermut, int *arows,
                             unsigned int nnzb, index2 *bi,
                             unsigned int *bpermut, int *bcols,
                             index2 *abi, int *abpos, index2 *aikbkj );

/* slower, but requiring less memory multiplication algorithm, */
/* the procedures with headers given below only find the distribution */
/* of nonzero product coefficients; their number must be determined */
/* prior to the call, using the appropriate *Count* procedure above */
boolean pkn_SPMmultMMCempty ( int nra, int nca, int ncb,
                              unsigned int nnza, index2 *ai,
                              unsigned int *apermut, int *acols, boolean ca,
                              unsigned int nnzb, index2 *bi,
                              unsigned int *bpermut, int *bcols, boolean cb,
                              index2 *abi );
boolean pkn_SPMmultMMTCempty ( int nra, int nca, int nrb,
                               unsigned int nnza, index2 *ai,
                               unsigned int *apermut, int *acols, boolean ca,
                               unsigned int nnzb, index2 *bi,
                               unsigned int *bpermut, int *brows, boolean rb,
                               index2 *abi );
boolean pkn_SPMmultMTMCempty ( int nra, int nca, int ncb,
                               unsigned int nnza, index2 *ai,
                               unsigned int *apermut, int *arows, boolean ra,
                               unsigned int nnzb, index2 *bi,
                               unsigned int *bpermut, int *bcols, boolean ba,
                               index2 *abi );

/* procedures for submatrices */
/* auxiliary procedures */
boolean pkn_SPsubMSortByRows ( unsigned int nrows, unsigned int ncols,
                               unsigned int nnz, index3 *ai, unsigned int *permut );
boolean pkn_SPsubMSortByCols ( unsigned int nrows, unsigned int ncols,
                               unsigned int nnz, index3 *ai, unsigned int *permut );
boolean pkn_SPsubMFindRows ( unsigned int nrows, unsigned int ncols, unsigned int nnz,
                             index3 *ai, unsigned int *permut, boolean ro,
                             int *rows );
boolean pkn_SPsubMFindCols ( unsigned int nrows, unsigned int ncols, unsigned int nnz,
                             index3 *ai, unsigned int *permut, boolean co,
                             int *cols );

/* fast multiplication algorithm - uses considerable additional memory */
boolean pkn_SPsubMCountMMnnzR ( int nra, int nca, int ncb,
                                unsigned int nnza, index3 *ai,
                                unsigned int *apermut, int *arows, boolean ra,
                                unsigned int nnzb, index3 *bi,
                                unsigned int *bpermut, int *brows, boolean rb,
                                unsigned int *nnzab, unsigned int *nmultab );
boolean pkn_SPsubMFindMMnnzR ( int nra, int nca, int ncb,
                               unsigned int nnza, index3 *ai,
                               unsigned int *apermut, int *arows,
                               unsigned int nnzb, index3 *bi,
                               unsigned int *bpermut, int *brows,
                               index2 *abi, int *abpos, index2 *aikbkj );
boolean pkn_SPsubMCountMMnnzC ( int nra, int nca, int ncb,
                                unsigned int nnza, index3 *ai,
                                unsigned int *apermut, int *acols, boolean ca,
                                unsigned int nnzb, index3 *bi,
                                unsigned int *bpermut, int *bcols, boolean cb,
                                unsigned int *nnzab, unsigned int *nmultab );
boolean pkn_SPsubMFindMMnnzC ( int nra, int nca, int ncb,
                               unsigned int nnza, index3 *ai,
                               unsigned int *apermut, int *acols,
                               unsigned int nnzb, index3 *bi,
                               unsigned int *bpermut, int *bcols,
                               index2 *abi, int *abpos, index2 *aikbkj );

boolean pkn_SPsubMCountMMTnnzR ( int nra, int nca, int nrb,
                                 unsigned int nnza, index3 *ai,
                                 unsigned int *apermut, int *arows, boolean ra,
                                 unsigned int nnzb, index3 *bi,
                                 unsigned int *bpermut, int *bcols, boolean cb,
                                 unsigned int *nnzab, unsigned int *nmultab );
boolean pkn_SPsubMFindMMTnnzR ( int nra, int nca, int nrb,
                                unsigned int nnza, index3 *ai,
                                unsigned int *apermut, int *arows,
                                unsigned int nnzb, index3 *bi,
                                unsigned int *bpermut, int *bcols,
                                index2 *abi, int *abpos, index2 *aikbkj );
boolean pkn_SPsubMCountMMTnnzC ( int nra, int nca, int nrb,
                                 unsigned int nnza, index3 *ai,
                                 unsigned int *apermut, int *acols, boolean ca,
                                 unsigned int nnzb, index3 *bi,
                                 unsigned int *bpermut, int *brows, boolean rb,
                                 unsigned int *nnzab, unsigned int *nmultab );
boolean pkn_SPsubMFindMMTnnzC ( int nra, int nca, int nrb,
                                unsigned int nnza, index3 *ai,
                                unsigned int *apermut, int *acols,
                                unsigned int nnzb, index3 *bi,
                                unsigned int *bpermut, int *brows,
                                index2 *abi, int *abpos, index2 *aikbkj );

boolean pkn_SPsubMCountMTMnnzR ( int nra, int nca, int ncb,
                                 unsigned int nnza, index3 *ai,
                                 unsigned int *apermut, int *acols, boolean ca,
                                 unsigned int nnzb, index3 *bi,
                                 unsigned int *bpermut, int *brows, boolean rb,
                                 unsigned int *nnzab, unsigned int *nmultab );
boolean pkn_SPsubMFindMTMnnzR ( int nra, int nca, int ncb,
                                unsigned int nnza, index3 *ai,
                                unsigned int *apermut, int *acols,
                                unsigned int nnzb, index3 *bi,
                                unsigned int *bpermut, int *brows,
                                index2 *abi, int *abpos, index2 *aikbkj );
boolean pkn_SPsubMCountMTMnnzC ( int nra, int nca, int ncb,
                                 unsigned int nnza, index3 *ai,
                                 unsigned int *apermut, int *arows, boolean ra,
                                 unsigned int nnzb, index3 *bi,
                                 unsigned int *bpermut, int *bcols, boolean cb,
                                 unsigned int *nnzab, unsigned int *nmultab );
boolean pkn_SPsubMFindMTMnnzC ( int nra, int nca, int ncb,
                                int unsigned nnza, index3 *ai,
                                unsigned int *apermut, int *arows,
                                unsigned int nnzb, index3 *bi,
                                unsigned int *bpermut, int *bcols,
                                index2 *abi, int *abpos, index2 *aikbkj );

#ifdef __cplusplus
}
#endif

#endif /* PKNUM_H*/

