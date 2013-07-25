
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* changes:                                                                  */
/* 22.07.2013, R. Putanowicz - static pointers to application reading        */
/*   procedures replaced by pointers in a bsf_UserReaders structure passed   */
/*   by a parameter to the reading procedure.                                */
/* 24.07.2013, P. Kiciak - changes related with integration of the above     */
/*   with the package.                                                       */

/* Header file for the libbsfile library - reading and writing text files */
/* with geometric data */

#ifndef BSFILE_H
#define BSFILE_H

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif
#ifndef PKNUM_H
#include "pknum.h"
#endif
#ifndef PKGEOM_H
#include "pkgeom.h"
#endif
#ifndef MULTIBS_H
#include "multibs.h"
#endif
#ifndef BSMESH_H
#include "bsmesh.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define BSF_MAX_NAME_LENGTH 64

#define BSF_SYMB_EOF       0
#define BSF_SYMB_ERROR     1
#define BSF_SYMB_INTEGER   2
#define BSF_SYMB_FLOAT     3
#define BSF_SYMB_LBRACE    4
#define BSF_SYMB_RBRACE    5
#define BSF_SYMB_PLUS      6
#define BSF_SYMB_MINUS     7
#define BSF_SYMB_STRING    8
#define BSF_SYMB_COMMA     9

#define BSF_FIRST_KEYWORD  10  /* number of the first keyword */
#define BSF_SYMB_BCURVE    10  /* consecutive numbers of keywords */
#define BSF_SYMB_BPATCH    11  /* sorted alphabetically */
#define BSF_SYMB_BSCURVE   12
#define BSF_SYMB_BSHOLE    13
#define BSF_SYMB_BSMESH    14
#define BSF_SYMB_BSPATCH   15
#define BSF_SYMB_CLOSED    16
#define BSF_SYMB_CPOINTS   17
#define BSF_SYMB_CPOINTSMK 18
#define BSF_SYMB_DEGREE    19
#define BSF_SYMB_DIM       20
#define BSF_SYMB_DOMAIN    21
#define BSF_SYMB_FACETS    22
#define BSF_SYMB_HALFEDGES 23
#define BSF_SYMB_KNOTS     24
#define BSF_SYMB_KNOTS_U   25
#define BSF_SYMB_KNOTS_V   26
#define BSF_SYMB_NAME      27
#define BSF_SYMB_RATIONAL  28
#define BSF_SYMB_SIDES     29
#define BSF_SYMB_UNIFORM   30
#define BSF_SYMB_VERTICES  31

#define BSF_NKEYWORDS 22
extern const char *bsf_keyword[BSF_NKEYWORDS];

extern FILE   *bsf_input, *bsf_output;

extern int    bsf_nextsymbol;
extern int    bsf_nextint;
extern double bsf_nextfloat;
extern char   *bsf_namebuffer;

/* ////////////////////////////////////////////////////////////////////////// */
boolean bsf_OpenInputFile ( const char *filename );
void bsf_CloseInputFile ( void );
void bsf_GetNextSymbol ( void );

void bsf_PrintErrorLocation ( void );

boolean bsf_ReadIntNumber ( int *number );
boolean bsf_ReadDoubleNumber ( double *number );
boolean bsf_ReadPointd ( int maxcpdimen, double *point, int *cpdimen );
int bsf_ReadPointsd ( int maxcpdimen, int maxnpoints,
                      double *points, int *cpdimen );
int bsf_ReadPointsMK ( int maxnpoints, byte *mk );

boolean bsf_ReadSpaceDim ( int maxdim, int *spdimen );
boolean bsf_ReadCurveDegree ( int maxdeg, int *degree );
boolean bsf_ReadPatchDegree ( int maxdeg, int *udeg, int *vdeg );

boolean bsf_ReadKnotSequenced ( int maxlastknot, int *lastknot, double *knots,
                                boolean *closed );

boolean bsf_ReadBezierCurve4d ( int maxdeg, int *deg, point4d *cpoints,
                                int *spdimen, boolean *rational,
                                byte *mk, char *name );

boolean bsf_ReadBSplineCurve4d ( int maxdeg, int maxlastknot, int maxncpoints,
                                 int *deg, int *lastknot, double *knots,
                                 boolean *closed, point4d *cpoints,
                                 int *spdimen, boolean *rational,
                                 byte *mk, char *name );

boolean bsf_ReadBezierPatch4d ( int maxdeg, int *udeg, int *vdeg,
                                int *pitch, point4d *cpoints,
                                int *spdimen, boolean *rational,
                                byte *mk, char *name );

boolean bsf_ReadBSplinePatch4d ( int maxdeg, int maxlastknot, int maxncpoints,
                                 int *udeg, int *lastknotu, double *knotsu,
                                 int *vdeg, int *lastknotv, double *knotsv,
                                 boolean *closed_u, boolean *closed_v,
                                 int *pitch, point4d *cpoints,
                                 int *spdimen, boolean *rational,
                                 byte *mk, char *name );

boolean bsf_ReadBSMesh4d ( int maxnv, int maxnhe, int maxnfac,
                           int *degree,
                           int *nv, BSMvertex *mv, int *mvhei, point4d *vc,
                           int *nhe, BSMhalfedge *mhe,
                           int *nfac, BSMfacet *mfac, int *mfhei,
                           int *spdimen, boolean *rational,
                           byte *mkv, char *name );

boolean bsf_ReadBSplineHoled ( int maxk, int *hole_k, double *knots,
                               point2d *domain_cp, point4d *hole_cp,
                               int *spdimen, boolean *rational,
                               byte *mk, char *name );

/* Below are defined 6 typedefs for function pointer. These typedef's
   make easier to hande user defined callbacks for reading data from
   BSTools files to user specified data structures.
   The typedefs are:
      bsf_BC_fptr  -- for Bezier Curves callback
      bsf_BSC_fptr -- for BSpline Curves callback
      bsf_BP_fptr  -- for Bezier Patch callback
      bsf_BSP_fptr -- for BSpline Patch callback
      bsf_BSM_fptr -- for BS Mesh callback
      bsf_BSH_fptr -- for BSpline Hole callback
*/
typedef void (*bsf_BC_fptr) ( void *userData, const char *name,
                              int degree, const point4d *cpoints,
                              int spdimen, boolean rational, byte *mk );

typedef void (*bsf_BSC_fptr) ( void *userData, const char *name,
                               int degree, int lastknot, const double *knots,
                               boolean closed,
                               const point4d *cpoints,
                               int spdimen, boolean rational, byte *mk );
typedef void (*bsf_BP_fptr) ( void *userData, const char *name, int udeg, int vdeg,
                              int pitch, const point4d *cpoints,
                              int spdimen, boolean rational, byte *mk );
typedef void (*bsf_BSP_fptr) ( void *userData, const char *name,
                               int udeg, int lastknotu, const double *knotsu,
                               int vdeg, int lastknotv, const double *knotsv,
                               boolean closed_u, boolean closed_v,
                               int pitch, const point4d *cpoints,
                               int spdimen, boolean rational, byte *mk );
typedef void (*bsf_BSM_fptr) ( void *userData, const char *name, int degree,
                               int nv, const BSMvertex *mv, const int *mvhei,
                               const point4d *vc,
                               int nhe, const BSMhalfedge *mhe,
                               int nfac, const BSMfacet *mfac, const int *mfhei,
                               int spdimen, boolean rational, byte *mkv );

typedef void (*bsf_BSH_fptr) ( void *userData, const char *name, int hole_k,
                               const double *knots,
                               const point2d *domain_cp,
                               const point4d *hole_cp,
                               int spdimen, boolean rational, byte *mk );

typedef struct {
  bsf_BC_fptr  BezierCurveReader;
  bsf_BSC_fptr BSplineCurveReader;
  bsf_BP_fptr  BezierPatchReader;
  bsf_BSP_fptr BSplinePatchReader;
  bsf_BSM_fptr BSMeshReader;
  bsf_BSH_fptr BSplineHoleReader;
  void *userData;   /* pointer to user data, bsfile library does not
                       use it, it is just passed to reader callbacks */
  int bc_maxdeg;       /* maximal degree of Bezier curves */
  int bsc_maxdeg;      /* maximal degree of B-spline curves */
  int bsc_maxlkn;      /* maximum number of the last knot for B-spline curves */
  int bp_maxdeg;       /* maximal degree of Bezier patches */
  int bsp_maxdeg;      /* maximal degree of B-spline patches */
  int bsp_maxlkn;      /* maximal number of the last knot for B-spline patches */
  int bsm_maxdeg;      /* maximal degree of mesh-represented surfaces */
  int bsm_maxnv;       /* maximum number of mesh vertices */
  int bsm_maxnhe;      /* maximum number of mesh halfedges */
  int bsm_maxnfac;     /* maximum number of mesh facets */ 
} bsf_UserReaders;

/* Allocate on stack new user readers data structure and initialize 
   it to be empty (i.e. all callback pointers and user data pointer
   are null, all data size restrictions are set some default values).
 */
bsf_UserReaders * bsf_NewReaders ( void );

/* Set up readers data structure to an empty state (see bsf_NewReaders)
 */
void bsf_ClearReaders ( bsf_UserReaders *readers );

#define bsf_DeleteReaders(readers) \
do { \
  PKV_FREE(readers); \
} while(0) 

void bsf_BC4ReadFuncd ( bsf_UserReaders *readers, bsf_BC_fptr BCReader,
                        int maxdeg );
void bsf_BSC4ReadFuncd ( bsf_UserReaders *readers, bsf_BSC_fptr BSCReader,
                         int maxdeg, int maxlastknot );
void bsf_BP4ReadFuncd ( bsf_UserReaders *readers, bsf_BP_fptr BPReader,
                        int maxdeg );
void bsf_BSP4ReadFuncd ( bsf_UserReaders *readers, bsf_BSP_fptr BSPReader,
                         int maxdeg, int maxlastknot );
void bsf_BSM4ReadFuncd ( bsf_UserReaders *readers, bsf_BSM_fptr BSMReader,
                         int maxdeg, int maxnv, int maxnhe, int maxnfac );
void bsf_BSH4ReadFuncd ( bsf_UserReaders *readers, bsf_BSH_fptr BSHReader );

boolean bsf_ReadBSFiled ( const char *filename, bsf_UserReaders *readers );

/* ////////////////////////////////////////////////////////////////////////// */
boolean bsf_OpenOutputFile ( char *filename, boolean append );
void bsf_CloseOutputFile ( void );

void bsf_WriteComment ( char *comment );
void bsf_WriteDoubleNumber ( double x );
void bsf_WritePointd ( int cpdimen, const double *point );
void bsf_WritePointsd ( int cpdimen, int cols, int rows, int pitch,
                        const double *points );
void bsf_WritePointsMK ( int npoints, const byte *mk );

void bsf_WriteSpaceDim ( int spdimen, boolean rational );
void bsf_WriteCurveDegree ( int degree );
void bsf_WritePatchDegree ( int udeg, int vdeg );

void bsf_WriteKnotSequenced ( int lastknot, const double *knots, boolean closed );

boolean bsf_WriteBezierCurved ( int spdimen, int cpdimen, boolean rational,
                                int deg, const double *cpoints,
                                const byte *mk, const char *name );

boolean bsf_WriteBSplineCurved ( int spdimen, int cpdimen, boolean rational,
                                 int deg, int lastknot, const double *knots,
                                 boolean closed,
                                 const double *cpoints, const byte *mk,
                                 const char *name );

boolean bsf_WriteBezierPatchd ( int spdimen, int cpdimen, boolean rational,
                                int udeg, int vdeg,
                                int pitch, const double *cpoints, const byte *mk,
                                const char *name );

boolean bsf_WriteBSplinePatchd ( int spdimen, int cpdimen, boolean rational,
                                 int udeg, int lastknotu, const double *knotsu,
                                 int vdeg, int lastknotv, const double *knotsv,
                                 boolean closed_u, boolean closed_v,
                                 int pitch, const double *cpoints, const byte *mk,
                                 const char *name );

boolean bsf_WriteBSMeshd ( int spdimen, int cpdimen, boolean rational, int degree,
                           int nv, const BSMvertex *mv, const int *mvhei,
                           const double *vc,
                           int nhe, const BSMhalfedge *mhe,
                           int nfac, const BSMfacet *mfac, const int *mfhei,
                           const byte *mkv, const char *name );

boolean bsf_WriteBSplineHoled ( int hole_k, const double *knots,
                                const point2d *domain_cp,
                                const point3d *hole_cp, const byte *mk,
                                const char *name );

#ifdef __cplusplus
}
#endif

#endif /*BSFILE_H*/

