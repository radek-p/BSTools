
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/times.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "bsmesh.h"
#include "bsfile.h"
#include "egholed.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"

/* ////////////////////////////////////////////////////////////////////////// */
/* some limits for the mesh size */
#define MAXNV    30000
#define MAXNHE  120000
#define MAXNFAC  30000

#define C    0.01
#define NKN1 6
#define NKN2 8
/* ////////////////////////////////////////////////////////////////////////// */
/* file names */
char *ifn = NULL;  /* input mesh  */
char *ofn = NULL;  /* output mesh */
char *cfn = NULL;  /* coarse mesh */

/* options */
int     nbl = 0;
boolean shape_only     = false;
boolean no_constraints = false;
boolean write_all      = false;
boolean verbose        = false;
boolean alternative    = false;

/* status */
boolean result_written = false;

/* representation of the mesh to optimise */
int         nv, nhe, nfac, *mvhei = NULL, *mfhei = NULL;
BSMvertex   *mv = NULL;
point3d     *mvcp = NULL;
BSMhalfedge *mhe = NULL;
BSMfacet    *mfac = NULL;
byte        *mkcp = NULL;
boolean     mesh_ok = false;

/* representation of the coarse mesh to define a~preconditioner */
int         cnv, cnhe, cnfac, *cmvhei = NULL, *cmfhei = NULL;
BSMvertex   *cmv = NULL;
BSMhalfedge *cmhe = NULL;
BSMfacet    *cmfac = NULL;
boolean     cmesh_ok = false;

/* ////////////////////////////////////////////////////////////////////////// */
/* parsing the command line and printing usage instruction */
boolean ReadCmdLineParameters ( int argc, char **argv )
{
  boolean is_input, is_output, is_coarse;
  int     i;

  is_input = is_output = is_coarse = false;
  for ( i = 1; i < argc; i++ ) {
    if ( !strcmp ( "-i", argv[i] ) ) {
      ifn = argv[++i];
      is_input = true;
    }
    else if ( !strcmp ( "-o", argv[i] ) ) {
      if ( is_output )
        goto failure;
      ofn = argv[++i];
      is_output = true;
    }
    else if ( !strcmp ( "-s", argv[i] ) ) {
      if ( nbl )
        goto failure;
      shape_only = true;
    }
    else if ( !strcmp ( "-c", argv[i] ) ) {
      if ( is_coarse )
        goto failure;
      cfn = argv[++i];
      is_coarse = true;
    }
    else if ( !strcmp ( "-b", argv[i] ) ) {
      if ( shape_only )
        goto failure;
      if ( sscanf ( argv[++i], "%d", &nbl ) != 1 )
        goto failure;
      if ( nbl < 2 || nbl > 32 )
        goto failure;
    }
    else if ( !strcmp ( "-q", argv[i] ) )
      alternative = true;
    else if ( !strcmp ( "-noconstr", argv[i] ) )
      no_constraints = true;
    else if ( !strcmp ( "-e", argv[i] ) )
      write_all = true;
    else if ( !strcmp ( "-v", argv[i] ) )
      verbose = true;
    else {
      printf ( "Invalid option: %s\n", argv[i] );
      goto failure;
    }
  }
  if ( !is_input || !is_output )
    goto failure;
  return true;

failure:
  return false;  
} /*ReadCmdLineParameters*/

void PrintUsageInfo ( void )
{
  printf ( "Usage:\n\n" );
  printf ( "  optblmesh [options] -i input_mesh.bs -o output_mesh.bs\n\n" );
  printf ( "Options:\n" );
  printf ( "  -a         (alternative multilevel algorithm)\n" );
  printf ( "  -b n       (two-level algorithm,\n" );
  printf ( "              n is the number of blocks, from 2 to 32)\n" );
  printf ( "  -c coarse_mesh.bs\n" );
  printf ( "  -e         (write the mesh after every iteration)\n" );
  printf ( "  -noconstr  (ignore the constraints)\n" );
  printf ( "  -s         (shape optimization, option -b must not be given)\n" );
  printf ( "  -v         (print more information)\n\n" );
} /*PrintUsageInfo*/

/* ////////////////////////////////////////////////////////////////////////// */
/* reading data files */
void StoreTheMesh ( const char *name, int degree,
                    int _nv, const BSMvertex *_mv, const int *_mvhei,
                    const point4d *_mvcp, int _nhe, const BSMhalfedge *_mhe,
                    int _nfac, const BSMfacet *_mfac, const int *_mfhei,
                    int spdimen, boolean rational, byte *_mkcp )
{
  int i;

        /* we take into account only the first mesh from the input file */
  if ( mesh_ok )
    return;
        /* the data must be copied from the scratch memory */
        /* to the arrays properly allocated */
  nv   = _nv;
  nhe  = _nhe;
  nfac = _nfac;
  PKV_MALLOC ( mv,    nv*sizeof(BSMvertex) );
  PKV_MALLOC ( mvhei, nhe*sizeof(int) );
  PKV_MALLOC ( mvcp,  nv*sizeof(point3d) );
  PKV_MALLOC ( mhe,   nhe*sizeof(BSMhalfedge) );
  PKV_MALLOC ( mfac,  nfac*sizeof(BSMfacet) );
  PKV_MALLOC ( mfhei, nhe*sizeof(int) );
  PKV_MALLOC ( mkcp,  nv*sizeof(byte) );
  if ( !mv || !mvhei || !mvcp || !mhe || !mfac || !mfhei || !mkcp ) {
    printf ( "Error: could not allocate arrays for the mesh\n" );
    exit ( 1 );
  }
  memcpy ( mv,    _mv,    nv*sizeof(BSMvertex) );
  memcpy ( mvhei, _mvhei, nhe*sizeof(int) );
  memcpy ( mhe,   _mhe,   nhe*sizeof(BSMhalfedge) );
  memcpy ( mfac,  _mfac,  nfac*sizeof(BSMfacet) );
  memcpy ( mfhei, _mfhei, nhe*sizeof(int) );
        /* Constraints may be imposed by marking vertices;       */
        /* the fixed vertices are the ones, which have nonzero   */
        /* bytes in the _mkcp array, which may be absent.        */
        /* In the ifn.bs file there should be a~section labelled */
        /* with the cpointsmk keyword, made of numbers. Nonfixed */
        /* vertices must have 0 in this section.                 */
  if ( _mkcp && !no_constraints ) {
    memcpy ( mkcp, _mkcp, nv*sizeof(byte) );
        /* as the program pozwalaj writes the merkings of vertices */
        /* with the least significant bit set, below we clear it.  */
        /* Thus fixed are the vertices with any other bit set.     */
    for ( i = 0; i < nv; i++ )
      mkcp[i] &= ~0x01;
  }
  else
    memset ( mkcp, 0, nv*sizeof(byte) );
  for ( i = 0; i < nv; i++ )
    Point4to3d ( &_mvcp[i], &mvcp[i] );
  mesh_ok = true;
} /*StoreTheMesh*/

boolean ReadInputMesh ( char *ifn )
{
        /* register the procedure to store the data read from the file */
  bsf_BSM4ReadFuncd ( StoreTheMesh, 10, MAXNV, MAXNHE, MAXNFAC );
  mesh_ok = false;
        /* read the file contents */
  bsf_ReadBSFiled ( ifn );
  if ( !mesh_ok )
    printf ( "Error: could not read a mesh from the input file.\n" );
  return mesh_ok;
} /*ReadInputMesh*/

void StoreTheCoarseMesh ( const char *name, int degree,
                    int _nv, const BSMvertex *_mv, const int *_mvhei,
                    const point4d *_mvcp, int _nhe, const BSMhalfedge *_mhe,
                    int _nfac, const BSMfacet *_mfac, const int *_mfhei,
                    int spdimen, boolean rational, byte *_mkcp )
{
        /* we take into account only the first mesh from the input file */
  if ( cmesh_ok )
    return;
        /* the data must be copied from the scratch memory */
        /* to the arrays properly allocated */
  cnv   = _nv;
  cnhe  = _nhe;
  cnfac = _nfac;
  PKV_MALLOC ( cmv,    cnv*sizeof(BSMvertex) );
  PKV_MALLOC ( cmvhei, cnhe*sizeof(int) );
  PKV_MALLOC ( cmhe,   cnhe*sizeof(BSMhalfedge) );
  PKV_MALLOC ( cmfac,  cnfac*sizeof(BSMfacet) );
  PKV_MALLOC ( cmfhei, cnhe*sizeof(int) );
  if ( !cmv || !cmvhei || !cmhe || !cmfac || !cmfhei ) {
    printf ( "Error: could not allocate arrays for the coarse mesh\n" );
    return;
  }
  memcpy ( cmv,    _mv,    cnv*sizeof(BSMvertex) );
  memcpy ( cmvhei, _mvhei, cnhe*sizeof(int) );
  memcpy ( cmhe,   _mhe,   cnhe*sizeof(BSMhalfedge) );
  memcpy ( cmfac,  _mfac,  cnfac*sizeof(BSMfacet) );
  memcpy ( cmfhei, _mfhei, cnhe*sizeof(int) );
  cmesh_ok = true;
} /*StoreTheCoarseMesh*/

void ReadCoarseMesh ( char *cfn )
{
  if ( !cfn )
    return;
        /* register the procedure to store the data read from the file */
  bsf_BSM4ReadFuncd ( StoreTheCoarseMesh, 10, MAXNV, MAXNHE, MAXNFAC );
  cmesh_ok = false;
        /* read the file contents */
  bsf_ReadBSFiled ( cfn );
        /* if reading the coarse mesh failed, it is still possible */
        /* to perform the optimisation, though the preconditioner  */
        /* used in this case is less effective */
  if ( !cmesh_ok )
    printf ( "Warning: could not read in the coarse mesh\n" );
} /*ReadCoarseMesh*/

/* ////////////////////////////////////////////////////////////////////////// */
/* writing out the result */
boolean WriteTheMesh ( char *ofn )
{
  if ( bsf_OpenOutputFile ( ofn, false ) ) {
    bsf_WriteBSMeshd ( 3, 3, false, 3, nv, mv, mvhei, &mvcp[0].x,
                       nhe, mhe, nfac, mfac, mfhei, mkcp, NULL );
    bsf_CloseOutputFile ();
    result_written = true;
    return true;
  }
  else
    return false;
} /*WriteTheMesh*/

/* ////////////////////////////////////////////////////////////////////////// */
/* optimisation */
boolean SetupRefinementMatrix ( int *m, int *n, int *nnz,
                                index2 **nnzi, double **nnzc )
{
  void        *sp, *sp1;
  int         nv1, nhe1, nfac1, *mvhei1, *mfhei1;
  BSMvertex   *mv1;
  BSMhalfedge *mhe1;
  BSMfacet    *mfac1;
  int         nnz1;
  index2      *nzi1;
  double      *nzc1;
  int         nv2, nhe2, nfac2, *mvhei2, *mfhei2;
  BSMvertex   *mv2;
  BSMhalfedge *mhe2;
  BSMfacet    *mfac2;
  int         nnz2;
  index2      *nzi2;
  double      *nzc2;
  int         nv0, nmult, *permut1, *cols1, *permut2, *cols2;

  sp = pkv_GetScratchMemTop ();

  mvhei1 = mfhei1 = mvhei2 = mfhei2 = NULL;
  mv1 = mfac1 = mv2 = mfac2 = NULL;
  mhe1 = mhe2 = NULL;
  nzi1 = nzi2 = NULL;
  nzc1 = nzc2 = NULL;
  
  if ( !mesh_ok || !cmesh_ok || nv <= cnv )
    goto failure;
        /* the first refinement of the coarse mesh */
  if ( !bsm_RefinementMatd ( 3, cnv, cmv, cmvhei, cnhe, cmhe, cnfac, cmfac, cmfhei,
                        &nv1, &mv1, &mvhei1, &nhe1, &mhe1, &nfac1, &mfac1, &mfhei1,
                        &nnz1, &nzi1, &nzc1 ) )
    goto failure;
        /* is it just one refinement step? if yes, we're done */
  if ( nv1 > nv )
    goto failure;
  if ( nv1 == nv ) {
    if ( nhe1 != nhe || nfac1 != nfac )
      goto failure;
    *nnz = nnz1;
    *nnzi = nzi1;
    *nnzc = nzc1;
    goto success;
  }
        /* no; it is then necessary to find the composition of a number */
        /* of refinement operations, i.e. the product of their matrices */
  nv0 = cnv;
  while ( nv1 < nv ) {
          /* get the next refinement operation */
    if ( !bsm_RefinementMatd ( 3, nv1, mv1, mvhei1, nhe1, mhe1, nfac1, mfac1, mfhei1,
                        &nv2, &mv2, &mvhei2, &nhe2, &mhe2, &nfac2, &mfac2, &mfhei2,
                        &nnz2, &nzi2, &nzc2 ) )
      goto failure;
    if ( nv2 > nv )
      goto failure;
          /* multiply the matrices in the sparse representations */
    sp1 = pkv_GetScratchMemTop ();
    permut1 = pkv_GetScratchMemi ( nnz1 );
    cols1   = pkv_GetScratchMemi ( nv0+1 );
    permut2 = pkv_GetScratchMemi ( nnz2 );
    cols2   = pkv_GetScratchMemi ( nv1+1 );
    if ( !permut1 || !cols1 || !permut2 || !cols2 )
      goto failure;
            /* counting the nonzero coefficients of the product and */
            /* mempry allocation */
    if ( !pkn_SPMCountMMnnzC ( nv2, nv1, nv0,
                               nnz2, nzi2, permut2, cols2, false,
                               nnz1, nzi1, permut1, cols1, false,
                               nnz, &nmult ) )
      goto failure;
    PKV_MALLOC ( *nnzi, (*nnz)*sizeof(index2) );
    PKV_MALLOC ( *nnzc, (*nnz)*sizeof(double) );
    if ( !(*nnzi) || !(*nnzc) )
      goto failure;
            /* actual multiplication */
    if ( !pkn_SPMmultMMCd ( nv2, nv1, nv0,
                            nnz2, nzi2, nzc2, permut2, cols2, true,
                            nnz1, nzi1, nzc1, permut1, cols1, true,
                            *nnzi, *nnzc ) )
      goto failure;
            /* we do not need the factors any more; deallocate them */
            /* as well as the coarser of the two meshes mesh */
    PKV_FREE ( mv1 );     mv1 = mv2;        mv2 = NULL;
    PKV_FREE ( mvhei1 );  mvhei1 = mvhei2;  mvhei2 = NULL;
    PKV_FREE ( mhe1 );    mhe1 = mhe2;      mhe2 = NULL;
    PKV_FREE ( mfac1 );   mfac1 = mfac2;    mfac2 = NULL;
    PKV_FREE ( mfhei1 );  mfhei1 = mfhei2;  mfhei2 = NULL;
    nv1 = nv2;
    nhe1 = nhe2;
    nfac1 = nfac2;
    PKV_FREE ( nzi1 );
    PKV_FREE ( nzc1 );
    PKV_FREE ( nzi2 );
    PKV_FREE ( nzc2 );
    pkv_SetScratchMemTop ( sp1 );
    if ( nv1 == nv ) {
      if ( nhe1 == nhe && nfac1 == nfac )
        goto success;
      else
        goto failure;
    }
    nnz1 = *nnz;
    nzi1 = *nnzi;
    nzc1 = *nnzc;
  }

success:
  if ( mv1 )    PKV_FREE ( mv1 );
  if ( mvhei1 ) PKV_FREE ( mvhei1 );
  if ( mhe1 )   PKV_FREE ( mhe1 );
  if ( mfac1 )  PKV_FREE ( mfac1 );
  if ( mfhei1 ) PKV_FREE ( mfhei1 );
  if ( mv2 )    PKV_FREE ( mv2 );
  if ( mvhei2 ) PKV_FREE ( mvhei2 );
  if ( mhe2 )   PKV_FREE ( mhe2 );
  if ( mfac2 )  PKV_FREE ( mfac2 );
  if ( mfhei2 ) PKV_FREE ( mfhei2 );
  if ( nzi1 )   PKV_FREE ( nzi1 );
  if ( nzc1 )   PKV_FREE ( nzc1 );
  if ( nzi2 )   PKV_FREE ( nzi2 );
  if ( nzc2 )   PKV_FREE ( nzc2 );
  *m = nv;
  *n = cnv;
  if ( mv1 )
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( mv1 )    PKV_FREE ( mv1 );
  if ( mvhei1 ) PKV_FREE ( mvhei1 );
  if ( mhe1 )   PKV_FREE ( mhe1 );
  if ( mfac1 )  PKV_FREE ( mfac1 );
  if ( mfhei1 ) PKV_FREE ( mfhei1 );
  if ( mv2 )    PKV_FREE ( mv2 );
  if ( mvhei2 ) PKV_FREE ( mvhei2 );
  if ( mhe2 )   PKV_FREE ( mhe2 );
  if ( mfac2 )  PKV_FREE ( mfac2 );
  if ( mfhei2 ) PKV_FREE ( mfhei2 );
  if ( nzi1 )   PKV_FREE ( nzi1 );
  if ( nzc1 )   PKV_FREE ( nzc1 );
  if ( nzi2 )   PKV_FREE ( nzi2 );
  if ( nzc2 )   PKV_FREE ( nzc2 );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*SetupRefinementMatrix*/

boolean OptimiseShapeOnly ( void )
{
  void    *d;
  int     minlev, maxlev;
  int     m, n, nnz;
  index2  *nzi;
  double  *nzc;
  boolean finished, result;

        /* preparation */
  nzi = NULL;
  nzc = NULL;
  if ( cmesh_ok ) {  /* try to use the coarse mesh preconditioner */
    if ( !g2mbl_MLCPSSuggestNLevels ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                      mkcp, &minlev, &maxlev ) )
      goto failure;
    if ( minlev < 2 )
      goto no_coarse;
    if ( !SetupRefinementMatrix ( &m, &n, &nnz, &nzi, &nzc ) )
      goto no_coarse;
    if ( !g2mbl_MLSCMPOptInitd ( nv, mv, mvhei, mvcp, nhe, mhe,
                                 nfac, mfac, mfhei, mkcp,
                                 cnv, nnz, nzi, nzc,
                                 NKN1, NKN2, minlev, &d ) )
      goto failure;
  }
  else {
no_coarse:
    if ( !g2mbl_MLSSuggestNLevels ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                    mkcp, &minlev, &maxlev ) )
      goto failure;
    if ( !g2mbl_MLSOptInitd ( nv, mv, mvhei, mvcp, nhe, mhe,
                              nfac, mfac, mfhei, mkcp,
                              NKN1, NKN2, minlev, &d ) )
      goto failure;
  }
  if ( verbose )
    g2mbl_MLSetLogLeveld ( d, 2 );
        /* actual optimisation */
  finished = false;
  do {
    if ( alternative )
      result = g2mbl_MLSCOptIterd ( d, &finished );
    else
      result = g2mbl_MLSOptIterd ( d, &finished );
    if ( !result ) {
      g2mbl_MLOptDeallocated ( &d );
      goto failure;
    }
    if ( write_all || finished ) {
      if ( !WriteTheMesh ( ofn ) )
        goto failure;
    }
  } while ( !finished );
  g2mbl_MLOptDeallocated ( &d );
  if ( nzi ) PKV_FREE ( nzi );
  if ( nzc ) PKV_FREE ( nzc );
  return true;

failure:
  if ( nzi ) PKV_FREE ( nzi );
  if ( nzc ) PKV_FREE ( nzc );
  return false;
} /*OptimiseShapeOnly*/

boolean OptimiseTwoLevel ( int nbl )
{
  void   *d;
  int    m, n, nnz;
  index2 *nzi;
  double *nzc;
  boolean finished;

        /* preparation */
  nzi = NULL;
  nzc = NULL;
  if ( cmesh_ok ) {
    if ( !SetupRefinementMatrix ( &m, &n, &nnz, &nzi, &nzc ) )
      goto no_coarse;
    if ( !g2mbl_InitBlCMPSurfaceOptd ( nv, mv, mvhei, mvcp, nhe, mhe,
                                       nfac, mfac, mfhei, mkcp,
                                       cnv, nnz, nzi, nzc,
                                       C, -1.0, -1.0, NKN1, NKN2, nbl, &d ) )
      goto failure;
  }
  else {
no_coarse:
    if ( !g2mbl_InitBlSurfaceOptAltBLMTd ( nv, mv, mvhei, mvcp, nhe, mhe,
                                           nfac, mfac, mfhei, mkcp,
                                           C, -1.0, -1.0, NKN1, NKN2, nbl, &d ) )
      goto failure;
  }
        /* actual optimisation */
  finished = false;
  do {
    if ( !g2mbl_IterBlSurfaceOptAltBLMTd ( d, &finished ) ) {
      g2mbl_OptLMTDeallocated ( &d );
      goto failure;
    }
    if ( write_all || finished ) {
      if ( !WriteTheMesh ( ofn ) )
        goto failure;
    }
  } while ( !finished );


  if ( nzi ) PKV_FREE ( nzi );
  if ( nzc ) PKV_FREE ( nzc );
  return true;

failure:
  if ( nzi ) PKV_FREE ( nzi );
  if ( nzc ) PKV_FREE ( nzc );
  return false;
} /*OptimiseTwoLevel*/

boolean OptimiseMultiLevel ( void )
{
  void    *d;
  int     minlev, maxlev;
  int     m, n, nnz;
  index2  *nzi;
  double  *nzc;
  boolean finished, result;

        /* preparation */
  nzi = NULL;
  nzc = NULL;
  if ( cmesh_ok ) {  /* try to use the coarse mesh preconditioner */
    if ( !g2mbl_MLCPSuggestNLevels ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                     mkcp, &minlev, &maxlev ) )
      goto failure;
    if ( minlev < 2 )
      goto no_coarse;
    if ( !SetupRefinementMatrix ( &m, &n, &nnz, &nzi, &nzc ) )
      goto no_coarse;
    if ( !g2mbl_MLCMPOptInitd ( nv, mv, mvhei, mvcp, nhe, mhe,
                                nfac, mfac, mfhei, mkcp,
                                cnv, nnz, nzi, nzc,
                                C, -1.0, -1.0, NKN1, NKN2, minlev, &d ) )
      goto failure;
  }
  else {
no_coarse:
    if ( !g2mbl_MLSuggestNLevels ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                   mkcp, &minlev, &maxlev ) )
      goto failure;
    if ( !g2mbl_MLOptInitd ( nv, mv, mvhei, mvcp, nhe, mhe,
                             nfac, mfac, mfhei, mkcp,
                             C, -1.0, -1.0, NKN1, NKN2, minlev, &d ) )
      goto failure;
  }
  if ( verbose )
    g2mbl_MLSetLogLeveld ( d, 2 );
        /* actual optimisation */
  finished = false;
  do {
    if ( alternative )
      result = g2mbl_MLSOptIterd ( d, &finished );
    else
      result = g2mbl_MLOptIterd ( d, &finished );
    if ( !result ) {
      g2mbl_MLOptDeallocated ( &d );
      goto failure;
    }
    if ( write_all || finished ) {
      if ( !WriteTheMesh ( ofn ) )
        goto failure;
    }
  } while ( !finished );

  g2mbl_MLOptDeallocated ( &d );
  if ( nzi ) PKV_FREE ( nzi );
  if ( nzc ) PKV_FREE ( nzc );
  return true;

failure:
  if ( nzi ) PKV_FREE ( nzi );
  if ( nzc ) PKV_FREE ( nzc );
  return false;
} /*OptimiseMultiLevel*/

boolean Optimise ( void )
{
  boolean    result;
  struct tms start, stop;

  times ( &start );
        /* choose the proper algorithm */
  if ( shape_only )
    result = OptimiseShapeOnly ();
  else if ( nbl )
    result = OptimiseTwoLevel ( nbl );
  else
    result = OptimiseMultiLevel ();
  times ( &stop );
  printf ( "total time = %6.2f\n",
           (double)(stop.tms_utime-start.tms_utime)/
           (double)sysconf(_SC_CLK_TCK) );
        /* After the computations there are data structures left,   */
        /* which describe the basis functions in the special domain */
        /* elements, and which might be reused (this would save     */
        /* a little time), but they should be deallocated, if the   */
        /* application has other things to do and needs the memory. */
        /* So, we deallocate, why not. */
  g2mbl_CleanupHoleDomainsd ();
  return result;
} /*Optimise*/

/* ////////////////////////////////////////////////////////////////////////// */
/* deallocate the memory occupied by the representation of the meshes */
void DeallocateMeshes ( void )
{
  if ( mvhei  ) PKV_FREE ( mvhei );
  if ( mfhei  ) PKV_FREE ( mfhei );
  if ( mv     ) PKV_FREE ( mv );
  if ( mvcp   ) PKV_FREE ( mvcp );
  if ( mhe    ) PKV_FREE ( mhe );
  if ( mfac   ) PKV_FREE ( mfac );
  if ( mkcp   ) PKV_FREE ( mkcp );
  if ( cmvhei ) PKV_FREE ( cmvhei );
  if ( cmfhei ) PKV_FREE ( cmfhei );
  if ( cmv    ) PKV_FREE ( cmv );
  if ( cmhe   ) PKV_FREE ( cmhe );
  if ( cmfac  ) PKV_FREE ( cmfac );
} /*DeallocateMeshes*/

/* ////////////////////////////////////////////////////////////////////////// */
/* main program - calls the procedures above */
int main ( int argc, char **argv )
{
  setvbuf ( stdout, NULL, _IONBF, 0 );
  pkv_InitScratchMem ( 256*1048576 ); /* allocate scratch memory pool, 256 MB */
  if ( !ReadCmdLineParameters ( argc, argv ) ) {
    PrintUsageInfo ();
    goto failure;
  }
  if ( !ReadInputMesh ( ifn ) )
    goto failure;
  ReadCoarseMesh ( cfn );
  if ( !Optimise () )
    goto failure;
  if ( result_written )
    printf ( "Result written in %s\n\n", ofn );
  else
    printf ( "Error: failed to write the result.\n\n" );
  DeallocateMeshes ();
  pkv_DestroyScratchMem ();
  exit ( 0 );

failure:
  DeallocateMeshes ();
  pkv_DestroyScratchMem ();
  exit ( 1 );
} /*main*/

