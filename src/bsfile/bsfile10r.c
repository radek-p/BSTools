
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* changes:                                                                  */
/* 22.07.2013, R. Putanowicz - static pointers to application reading        */
/*   procedures replaced by pointers in a bsf_UserReaders structure passed   */
/*   by a parameter to the reading procedure.                                */
/* 24.07.2013, P. Kiciak - changes related with integration of the above     */
/*   with the package.                                                       */
/* 27.12.2013, P. Kiciak - adding procedures reading colours and cameras.    */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "egholed.h"
#include "bsmesh.h"
#include "bsfile.h"

#include "bsfprivate.h"

/* ////////////////////////////////////////////////////////////////////////// */
static boolean _bsf_ReadBCurve ( bsf_UserReaders *readers )
{
  void    *sp;
  char    *name;
  point4d *cp;
  int     degree, spdimen;
  int     ident;
  boolean rational;

  sp = pkv_GetScratchMemTop ();
  name = pkv_GetScratchMem ( BSF_MAX_NAME_LENGTH+1 );
  cp = pkv_GetScratchMem ( (readers->bc_maxdeg+1)*sizeof(point4d) );
  if ( !name || !cp )
    goto failure;
  ident = -1;
  if ( !bsf_ReadBezierCurve4d ( readers->bc_maxdeg, &degree, cp,
                                &spdimen, &rational, name, &ident, readers ) )
    goto failure;
  if ( readers ) {
    if ( readers->BezierCurveReader )
      readers->BezierCurveReader ( readers->userData, name, ident, degree, cp,
                                   spdimen, rational );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_bsf_ReadBCurve*/

static boolean _bsf_ReadBPatch ( bsf_UserReaders *readers )
{
  void    *sp;
  char    *name;
  point4d *cp;
  int     udeg, vdeg, pitch, spdimen, maxdeg;
  int     ident;
  boolean rational;

  sp = pkv_GetScratchMemTop ();
  maxdeg = readers->bp_maxdeg;
  name = pkv_GetScratchMem ( BSF_MAX_NAME_LENGTH+1 );
  cp = pkv_GetScratchMem ( (maxdeg+1)*(maxdeg+1)*sizeof(point4d) );
  if ( !name || !cp )
    goto failure;
   ident = -1;
  if ( !bsf_ReadBezierPatch4d ( maxdeg, &udeg, &vdeg, &pitch, cp,
                                &spdimen, &rational, name, &ident, readers ) )
    goto failure;
  if ( readers ) {
    if ( readers->BezierPatchReader )
      readers->BezierPatchReader ( readers->userData, name, ident, udeg, vdeg, 
                                   pitch, cp, spdimen, rational );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_bsf_ReadBPatch*/

static boolean _bsf_ReadBSCurve ( bsf_UserReaders *readers )
{
  void    *sp;
  char    *name;
  point4d *cp;
  double  *knots;
  int     degree, lastknot, spdimen, maxdeg, maxlkn;
  int     ident;
  boolean rational, closed;

  sp = pkv_GetScratchMemTop ();
  maxlkn = readers->bsc_maxlkn;
  maxdeg = readers->bsc_maxdeg;
  name = pkv_GetScratchMem ( BSF_MAX_NAME_LENGTH+1 );
  cp = pkv_GetScratchMem ( maxlkn*sizeof(point4d) );
  knots = pkv_GetScratchMemd ( maxlkn+1 );
  if ( !name || !cp || !knots )
    goto failure;
   ident = -1;
  if ( !bsf_ReadBSplineCurve4d ( maxdeg, maxlkn, maxlkn,
                                 &degree, &lastknot, knots, &closed, cp,
                                 &spdimen, &rational, name, &ident, readers ) )
    goto failure;
  if ( readers ) {
    if ( readers->BSplineCurveReader )
      readers->BSplineCurveReader ( readers->userData, name, ident, degree, lastknot,
                                    knots, closed, cp, spdimen, rational );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_bsf_ReadBSCurve*/

static boolean _bsf_ReadBSPatch ( bsf_UserReaders *readers )
{
  void    *sp;
  char    *name;
  point4d *cp;
  double  *knotsu, *knotsv;
  int     spdimen, udeg, vdeg, lastknotu, lastknotv, pitch, maxdeg, maxlkn;
  int     ident;
  boolean closed_u, closed_v, rational;

  sp = pkv_GetScratchMemTop ();
  maxlkn = readers->bsp_maxlkn;
  maxdeg = readers->bsp_maxdeg;
  name = pkv_GetScratchMem ( BSF_MAX_NAME_LENGTH+1 );
  cp = pkv_GetScratchMem ( maxlkn*maxlkn*sizeof(point4d) );
  knotsu = pkv_GetScratchMemd ( maxlkn );
  knotsv = pkv_GetScratchMemd ( maxlkn );
  if ( !name || !cp || !knotsu || !knotsv )
    goto failure;
  ident = -1;
  if ( !bsf_ReadBSplinePatch4d ( maxdeg, maxlkn, maxlkn*maxlkn,
                                 &udeg, &lastknotu, knotsu,
                                 &vdeg, &lastknotv, knotsv,
                                 &closed_u, &closed_v, &pitch, cp,
                                 &spdimen, &rational, name, &ident, readers ) )
    goto failure;
  if ( readers ) {
    if ( readers->BSplinePatchReader )
      readers->BSplinePatchReader ( readers->userData,
               name, ident, udeg, lastknotu, knotsu, vdeg, lastknotv, knotsv,
               closed_u, closed_v, pitch, cp, spdimen, rational );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_bsf_ReadBSPatch*/

static boolean _bsf_ReadBSMesh ( bsf_UserReaders *readers )
{
  void        *sp;
  char        *name;
  BSMvertex   *mv;
  BSMhalfedge *mhe;
  BSMfacet    *mfac;
  int         *mvhei, *mfhei;
  point4d     *vc;
  int         nv, nhe, nfac, spdimen, degree, maxnv, maxnhe, maxnfac;
  int         ident;
  boolean     rational;

  sp = pkv_GetScratchMemTop ();
  maxnv = readers->bsm_maxnv;
  maxnhe = readers->bsm_maxnhe;
  maxnfac = readers->bsm_maxnfac;
  name = pkv_GetScratchMem ( BSF_MAX_NAME_LENGTH+1 );
  mv = pkv_GetScratchMem ( maxnv*sizeof(BSMvertex) );
  mhe = pkv_GetScratchMem ( maxnhe*sizeof(BSMhalfedge) );
  mfac = pkv_GetScratchMem ( maxnfac*sizeof(BSMfacet) );
  vc = pkv_GetScratchMem ( maxnv*sizeof(point4d) );
  mvhei = pkv_GetScratchMemi ( maxnhe );
  mfhei = pkv_GetScratchMemi ( maxnhe );
  if ( !name || !mv || !mhe || !mfac || !vc || !mvhei || !mfhei )
    goto failure;
  ident = -1;
  if ( !bsf_ReadBSMesh4d ( maxnv, maxnhe, maxnfac,
                           &degree, &nv, mv, mvhei, vc,
                           &nhe, mhe, &nfac, mfac, mfhei, &spdimen, &rational,
                           name, &ident, readers ) )
    goto failure;
  if ( readers ) {
    if ( readers->BSMeshReader )
      readers->BSMeshReader ( readers->userData, name, ident,
                              degree, nv, mv, mvhei, 
                              vc, nhe, mhe, nfac, mfac, mfhei,
                              spdimen, rational );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_bsf_ReadBSMesh*/

static boolean _bsf_ReadBSHole ( bsf_UserReaders *readers )
{
  void    *sp;
  int     hole_k, spdimen;
  char    *name;
  double  *knots;
  point2d *domcp;
  point4d *holecp;
  int     ident;
  boolean rational;

  sp = pkv_GetScratchMemTop ();
  name = pkv_GetScratchMem ( BSF_MAX_NAME_LENGTH+1 );
  knots = pkv_GetScratchMemd ( 11*GH_MAX_K );
  domcp = pkv_GetScratchMem ( (12*GH_MAX_K+1)*sizeof(point2d) );
  holecp = pkv_GetScratchMem ( (12*GH_MAX_K+1)*sizeof(point4d) );
  if ( !knots || !domcp || !holecp )
    goto failure;
  ident = -1;
  if ( !bsf_ReadBSplineHoled ( GH_MAX_K, &hole_k, knots, domcp, holecp,
                               &spdimen, &rational, name, &ident, readers ) )
    goto failure;
  if ( readers ) {
    if ( readers->BSplineHoleReader )
      readers->BSplineHoleReader ( readers->userData, name, ident,
                                   hole_k, knots, 
                                   domcp, holecp, spdimen, rational );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_bsf_ReadBSHole*/

static boolean _bsf_ReadPolyline ( bsf_UserReaders *readers )
{
  void    *sp;
  int     maxvert, nvert, spdimen;
  char    *name;
  point4d *vert;
  int     ident;
  boolean closed, rational;

  sp = pkv_GetScratchMemTop ();
  name = pkv_GetScratchMem ( BSF_MAX_NAME_LENGTH+1 );
  maxvert = readers->poly_maxvert;
  vert = pkv_GetScratchMem ( maxvert*sizeof(point4d) );
  if ( !name || !vert )
    goto failure;
  ident = -1;
  if ( !bsf_ReadPolyline4d ( maxvert, &nvert, vert, &spdimen,
                             &rational, &closed, name, &ident, readers ) )
    goto failure;
  if ( readers ) {
    if ( readers->PolylineReader )
      readers->PolylineReader ( readers->userData, name, ident,
                                nvert, vert, spdimen, closed, rational );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_bsf_ReadPolyline*/

boolean _bsf_ReadCPMark ( bsf_UserReaders *readers, int maxnpoints )
{
  void         *sp;
  unsigned int *mk;
  int          nmk;

  sp = pkv_GetScratchMemTop ();
  mk = (unsigned int*)pkv_GetScratchMemi ( maxnpoints );
  if ( !mk )
    goto failure;
  memset ( mk, 0, maxnpoints*sizeof(int) );
  nmk = bsf_ReadPointsMK ( maxnpoints, mk );
  if ( nmk <= 0 )
    goto failure;
  if ( readers ) {
    if ( readers->CPMarkReader )
      readers->CPMarkReader ( readers->userData, nmk, mk );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_bsf_ReadCPMark*/

boolean _bsf_ReadHEdgeMark ( bsf_UserReaders *readers, int maxnhe )
{
  void         *sp;
  unsigned int *mk;
  int          nmk;

  sp = pkv_GetScratchMemTop ();
  mk = (unsigned int*)pkv_GetScratchMemi ( maxnhe );
  if ( !mk )
    goto failure;
  memset ( mk, 0, maxnhe*sizeof(int) );
  nmk = bsf_ReadHEdgeMK ( maxnhe, mk );
  if ( nmk <= 0 )
    goto failure;
  if ( readers ) {
    if ( readers->HEdgeMarkReader )
      readers->HEdgeMarkReader ( readers->userData, nmk, mk );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_bsf_ReadHEdgeMark*/

boolean _bsf_ReadFacetMark ( bsf_UserReaders *readers, int maxnfac )
{
  void         *sp;
  unsigned int *mk;
  int          nmk;

  sp = pkv_GetScratchMemTop ();
  mk = (unsigned int*)pkv_GetScratchMemi ( maxnfac );
  if ( !mk )
    goto failure;
  memset ( mk, 0, maxnfac*sizeof(int) );
  nmk = bsf_ReadFacetMK ( maxnfac, mk );
  if ( nmk <= 0 )
    goto failure;
  if ( readers ) {
    if ( readers->FacetMarkReader )
      readers->FacetMarkReader ( readers->userData, nmk, mk );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_bsf_ReadFacetMark*/

boolean _bsf_ReadColour ( bsf_UserReaders *readers )
{
  point3d colour;

  if ( !bsf_ReadColour ( &colour ) )
    return false;
  if ( readers ) {
    if ( readers->ColourReader )
      readers->ColourReader ( readers->userData, &colour );
  }
  return true;
} /*_bsf_ReadColour*/

static boolean _bsf_ReadCamera ( bsf_UserReaders *readers )
{
  CameraRecd Camera;
  int        ident;

  ident = -1;
  if ( !bsf_ReadCamera ( &Camera, &ident ) )
    return false;
  if ( readers ) {
    if ( readers->CameraReader )
      readers->CameraReader ( readers->userData, ident, &Camera );
  }
  return true;
} /*_bsf_ReadCamera*/

/* procedure of reading an entire file; for each item the application */
/* routine is called in order to enter the data into the application */
/* data structures */
boolean bsf_ReadBSFiled ( const char *filename, bsf_UserReaders *readers )
{
  void    *sp;
  boolean signal_end;
  int     obj_type;

  sp = pkv_GetScratchMemTop ();
  if ( !bsf_OpenInputFile ( filename ) )
    return false;
  readers->done = false;
  for (;;) {
    obj_type = BSF_NONE;
    signal_end = false;
    if ( readers ) {
      if ( readers->BeginReader ) {
        switch ( bsf_nextsymbol ) {
    case BSF_SYMB_BCURVE:
          obj_type = BSF_BEZIER_CURVE;
          goto notify;
    case BSF_SYMB_BPATCH:
          obj_type = BSF_BEZIER_PATCH;
          goto notify;
    case BSF_SYMB_BSCURVE:
          obj_type = BSF_BSPLINE_CURVE;
          goto notify;
    case BSF_SYMB_BSPATCH:
          obj_type = BSF_BSPLINE_PATCH;
          goto notify;
    case BSF_SYMB_BSMESH:
          obj_type = BSF_BSPLINE_MESH;
          goto notify;
    case BSF_SYMB_BSHOLE:
          obj_type = BSF_BSPLINE_HOLE;
          goto notify;
    case BSF_SYMB_POLYLINE:
          obj_type = BSF_POLYLINE;
          goto notify;
    case BSF_SYMB_CAMERA:
          obj_type = BSF_CAMERA;
notify:
          readers->BeginReader ( readers->userData, obj_type );
          signal_end = true;
          break;
    default:
          break;
        }
      }
    }

    switch ( bsf_nextsymbol ) {
case BSF_SYMB_BCURVE:
      if ( !_bsf_ReadBCurve ( readers ) )
        goto failure;
      break;

case BSF_SYMB_BPATCH:
      if ( !_bsf_ReadBPatch ( readers ) )
        goto failure;
      break;

case BSF_SYMB_BSCURVE:
      if ( !_bsf_ReadBSCurve ( readers ) )
        goto failure;
      break;

case BSF_SYMB_BSPATCH:
      if ( !_bsf_ReadBSPatch ( readers ) )
        goto failure;
      break;

case BSF_SYMB_BSMESH:
      if ( !_bsf_ReadBSMesh ( readers ) )
        goto failure;
      break;

case BSF_SYMB_BSHOLE:
      if ( !_bsf_ReadBSHole ( readers ) )
        goto failure;
      break;

case BSF_SYMB_POLYLINE:
      if ( !_bsf_ReadPolyline ( readers ) )
        goto failure;
      break;

case BSF_SYMB_COLOR:
case BSF_SYMB_COLOUR:
      if ( !_bsf_ReadColour ( readers ) )
        goto failure;
      break;

case BSF_SYMB_CAMERA:
      if ( !_bsf_ReadCamera ( readers ) )
        goto failure;
      break;

case BSF_SYMB_EOF:
      goto finish;

default:
      goto failure;
    }
    if ( signal_end && readers->EndReader )
      readers->EndReader ( readers->userData, obj_type, true );
    if ( readers->done )
      goto finish;
  }

finish:
  bsf_CloseInputFile ();
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( signal_end && readers->EndReader )
    readers->EndReader ( readers->userData, obj_type, false );
  bsf_CloseInputFile ();
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsf_ReadBSFiled*/

void bsf_ClearReaders ( bsf_UserReaders *readers )
{
  if ( readers ) {
    readers->BeginReader        = NULL;
    readers->EndReader          = NULL;
    readers->BezierCurveReader  = NULL;
    readers->BSplineCurveReader = NULL;
    readers->BezierPatchReader  = NULL;
    readers->BSplinePatchReader = NULL;
    readers->BSMeshReader       = NULL;
    readers->BSplineHoleReader  = NULL;
    readers->PolylineReader     = NULL;
    readers->DepReader          = NULL;
    readers->TrimmedReader      = NULL;
    readers->CPMarkReader       = NULL;
    readers->HEdgeMarkReader    = NULL;
    readers->FacetMarkReader    = NULL;
    readers->CameraReader       = NULL;
    readers->ColourReader       = NULL;
    readers->userData           = NULL;
    readers->done = false;
        /* arbitrary default values */
    readers->bc_maxdeg    = 10;
    readers->bsc_maxdeg   = 10;
    readers->bsc_maxlkn   = 100;
    readers->bp_maxdeg    = 10;
    readers->bsp_maxdeg   = 10;
    readers->bsp_maxlkn   = 100;
    readers->bsm_maxdeg   = 10;
    readers->bsm_maxnv    = 1000;
    readers->bsm_maxnhe   = 2000;
    readers->bsm_maxnfac  = 1000;
    readers->poly_maxvert = 10;
    readers->maxdep       = 10;
  }
} /*bsf_ClearReaders*/

void bsf_BeginReadingFuncd ( bsf_UserReaders *readers,
                             bsf_BeginRead_fptr BeginReader )
{
  readers->BeginReader = BeginReader;
} /*bsf_BeginReadingFuncd*/

void bsf_EndReadingFuncd ( bsf_UserReaders *readers,
                           bsf_EndRead_fptr EndReader )
{
  readers->EndReader = EndReader;
} /*bsf_EndReadingFuncd*/

void bsf_BC4ReadFuncd ( bsf_UserReaders *readers, bsf_BC_fptr BCReader,
                        int maxdeg )
{
  readers->BezierCurveReader = BCReader;
  readers->bc_maxdeg = maxdeg;
} /*bsf_BC4ReadFuncd*/

void bsf_BSC4ReadFuncd ( bsf_UserReaders *readers, bsf_BSC_fptr BSCReader,
                         int maxdeg, int maxlastknot )
{
  readers->BSplineCurveReader = BSCReader;
  readers->bsc_maxdeg = maxdeg;
  readers->bsc_maxlkn = maxlastknot;
} /*bsf_BSC4ReadFuncd*/

void bsf_BP4ReadFuncd ( bsf_UserReaders *readers, bsf_BP_fptr BPReader,
                        int maxdeg )
{
  readers->BezierPatchReader = BPReader;
  readers->bp_maxdeg = maxdeg;
} /*bsf_BP4ReadFuncd*/

void bsf_BSP4ReadFuncd ( bsf_UserReaders *readers, bsf_BSP_fptr BSPReader,
                         int maxdeg, int maxlastknot )
{
  readers->BSplinePatchReader = BSPReader;
  readers->bsp_maxdeg = maxdeg;
  readers->bsp_maxlkn = maxlastknot;
} /*bsf_BSP4ReadFuncd*/

void bsf_BSM4ReadFuncd ( bsf_UserReaders *readers, bsf_BSM_fptr BSMReader,
                         int maxdeg, int maxnv, int maxnhe, int maxnfac )
{
  readers->BSMeshReader = BSMReader;
  readers->bsm_maxdeg = maxdeg;
  readers->bsm_maxnv = maxnv;
  readers->bsm_maxnhe = maxnhe;
  readers->bsm_maxnfac = maxnfac;
} /*bsf_BSM4ReadFuncd*/

void bsf_BSH4ReadFuncd ( bsf_UserReaders *readers, bsf_BSH_fptr BSHReader )
{
  readers->BSplineHoleReader = BSHReader;
} /*bsf_BSH4ReadFuncd*/

void bsf_Polyline4ReadFuncd ( bsf_UserReaders *readers, bsf_polyline_fptr PReader,
                              int maxvert )
{
  readers->PolylineReader = PReader;
  readers->poly_maxvert = maxvert;
} /*bsf_Polyline4ReadFuncd*/

void bsf_DependencyReadFunc ( bsf_UserReaders *readers,
                              bsf_dependency_fptr DepReader, int maxdep )
{
  readers->DepReader = DepReader;
  readers->maxdep = maxdep;
} /*bsf_DependencyReadFunc*/

void bsf_TrimmedReadFuncd ( bsf_UserReaders *readers,
                            bsf_trimmed_fptr TrimmedReader )
{
  readers->TrimmedReader = TrimmedReader;
} /*bsf_TrimmedReadFuncd*/

void bsf_CPMarkReadFunc ( bsf_UserReaders *readers,
                          bsf_CPMark_fptr CPMarkReader )
{
  readers->CPMarkReader = CPMarkReader;
} /*bsf_CPMarkReadFunc*/

void bsf_HalfedgeMarkReadFunc ( bsf_UserReaders *readers,
                                bsf_HEMark_fptr HEdgeMarkReader )
{
  readers->HEdgeMarkReader = HEdgeMarkReader;
} /*bsf_HEdgeMarkReadFunc*/

void bsf_FacetMarkReadFunc ( bsf_UserReaders *readers,
                             bsf_FacetMark_fptr FacetMarkReader )
{
  readers->FacetMarkReader = FacetMarkReader;
} /*bsf_FacetMarkReadFunc*/

void bsf_CameraReadFuncd ( bsf_UserReaders *readers,
                           bsf_Camera_fptr CameraReader )
{
  readers->CameraReader = CameraReader;
} /*bsf_CameraReadFuncd*/

void bsf_ColourReadFuncd ( bsf_UserReaders *readers,
                           bsf_Colour_fptr ColourReader )
{
  readers->ColourReader = ColourReader;
} /*bsf_ColourReadFuncd*/

