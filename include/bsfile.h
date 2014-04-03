
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* changes:                                                                  */
/* 22.07.2013, R. Putanowicz - static pointers to application reading        */
/*   procedures replaced by pointers in a bsf_UserReaders structure passed   */
/*   by a parameter to the reading procedure.                                */
/* 24.07.2013, P. Kiciak - changes related with integration of the above     */
/*   with the package.                                                       */
/* 7.01.2014, P. Kiciak - changes making it possible to extend the syntax    */
/*   of data files, to store object various attributes etc.                  */

/* Header file for the libbsfile library - reading and writing text files */
/* with geometric data */

#ifndef BSFILE_H
#define BSFILE_H

#include <stdio.h>

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
#ifndef CAMERA_H
#include "camera.h"
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

/* the subsequent keyword identifiers below must */
/* match the strings in the table in bsfile00r.c */
#define BSF_FIRST_KEYWORD    10  /* number of the first keyword */
#define BSF_SYMB_BCURVE      10  /* consecutive numbers of keywords */
#define BSF_SYMB_BPATCH      11  /* sorted alphabetically, with uppercase */
#define BSF_SYMB_BSCURVE     12  /* preceding all lowercase letters (ASCII) */
#define BSF_SYMB_BSHOLE      13
#define BSF_SYMB_BSMESH      14
#define BSF_SYMB_BSPATCH     15
#define BSF_SYMB_EULERANGLES 16
#define BSF_SYMB_CAMERA      17
#define BSF_SYMB_CLOSED      18
#define BSF_SYMB_COLOR       19  /* it can read both, AE and BE, */
#define BSF_SYMB_COLOUR      20  /* but it will write in BE */
#define BSF_SYMB_CPOINTS     21
#define BSF_SYMB_CPOINTSMK   22
#define BSF_SYMB_DEGREE      23
#define BSF_SYMB_DEPTH       24
#define BSF_SYMB_DIM         25
#define BSF_SYMB_DOMAIN      26
#define BSF_SYMB_FACETS      27
#define BSF_SYMB_FRAME       28
#define BSF_SYMB_HALFEDGES   29
#define BSF_SYMB_IDENT       30
#define BSF_SYMB_KNOTS       31
#define BSF_SYMB_KNOTS_U     32
#define BSF_SYMB_KNOTS_V     33
#define BSF_SYMB_NAME        34
#define BSF_SYMB_PARALLEL    35
#define BSF_SYMB_PERSPECTIVE 36
#define BSF_SYMB_POSITION    37
#define BSF_SYMB_RATIONAL    38
#define BSF_SYMB_SIDES       39
#define BSF_SYMB_UNIFORM     40
#define BSF_SYMB_VERTICES    41

/* namely, this table */
#define BSF_NKEYWORDS 32
extern const char *bsf_keyword[BSF_NKEYWORDS];

/* ////////////////////////////////////////////////////////////////////////// */
/* static variables */
extern FILE   *bsf_input, *bsf_output;

extern int    bsf_nextsymbol;
extern int    bsf_nextint;
extern double bsf_nextfloat;
extern char   *bsf_namebuffer;

extern int    bsf_current_indentation;

/* ////////////////////////////////////////////////////////////////////////// */
/* high-level reading procedures, read an entire file and ignore anything */
/* that an application is not interested in */

/* Below are defined 11 typedefs for function pointers. These typedef's   */
/* make easier to handle user defined callbacks for reading data from     */
/* BSTools files to user specified data structures.                       */
/* The typedefs are:                                                      */
/*    bsf_BeginRead_fptr -- called at the beginning of reading any object */
/*    bsf_EndRead_fptr   -- called at the end of reading an object        */
/*    bsf_BC_fptr        -- for Bezier Curves callback                    */
/*    bsf_BSC_fptr       -- for BSpline Curves callback                   */
/*    bsf_BP_fptr        -- for Bezier Patch callback                     */
/*    bsf_BSP_fptr       -- for BSpline Patch callback                    */
/*    bsf_BSM_fptr       -- for BS Mesh callback                          */
/*    bsf_BSH_fptr       -- for BSpline Hole callback                     */
/*    bsf_CPMark_fptr    -- for arrays of markings of control points      */
/*    bsf_Camera_fptr    -- for Camera object callback                    */
/*    bsf_Colour_fptr    -- for colour attribute callback                 */

typedef void (*bsf_BeginRead_fptr) ( void *userData, int obj_type );
typedef void (*bsf_EndRead_fptr) ( void *userData, boolean success );

typedef void (*bsf_BC_fptr) ( void *userData, const char *name, int ident,
                              int degree, const point4d *cpoints,
                              int spdimen, boolean rational );
typedef void (*bsf_BSC_fptr) ( void *userData, const char *name, int ident,
                               int degree, int lastknot, const double *knots,
                               boolean closed,
                               const point4d *cpoints,
                               int spdimen, boolean rational );
typedef void (*bsf_BP_fptr) ( void *userData, const char *name, int ident,
                              int udeg, int vdeg,
                              int pitch, const point4d *cpoints,
                              int spdimen, boolean rational );
typedef void (*bsf_BSP_fptr) ( void *userData, const char *name, int ident,
                               int udeg, int lastknotu, const double *knotsu,
                               int vdeg, int lastknotv, const double *knotsv,
                               boolean closed_u, boolean closed_v,
                               int pitch, const point4d *cpoints,
                               int spdimen, boolean rational );
typedef void (*bsf_BSM_fptr) ( void *userData, const char *name, int ident,
                               int degree,
                               int nv, const BSMvertex *mv, const int *mvhei,
                               const point4d *vc,
                               int nhe, const BSMhalfedge *mhe,
                               int nfac, const BSMfacet *mfac, const int *mfhei,
                               int spdimen, boolean rational );
typedef void (*bsf_BSH_fptr) ( void *userData, const char *name, int ident,
                               int hole_k, const double *knots,
                               const point2d *domain_cp,
                               const point4d *hole_cp,
                               int spdimen, boolean rational );
typedef void (*bsf_CPMark_fptr) ( void *userData, int ncp, unsigned int *mk );
typedef void (*bsf_Camera_fptr) ( void *userData, int ident, CameraRecd *Camera );
typedef void (*bsf_Colour_fptr) ( void *userData, point3d *colour );

typedef struct {
  void    *userData;   /* pointer to user data, bsfile library does not */
                       /* use it, it is just passed to reader callbacks */
  boolean done;        /* may be assigned true to stop reading */

                       /* data size limits */
  int     bc_maxdeg;   /* maximal degree of Bezier curves */
  int     bsc_maxdeg;  /* maximal degree of B-spline curves */
  int     bsc_maxlkn;  /* maximum number of the last knot for B-spline curves */
  int     bp_maxdeg;   /* maximal degree of Bezier patches */
  int     bsp_maxdeg;  /* maximal degree of B-spline patches */
  int     bsp_maxlkn;  /* maximal number of the last knot for B-spline patches */
  int     bsm_maxdeg;  /* maximal degree of mesh-represented surfaces */
  int     bsm_maxnv;   /* maximum number of mesh vertices */
  int     bsm_maxnhe;  /* maximum number of mesh halfedges */
  int     bsm_maxnfac; /* maximum number of mesh facets */ 

                       /* pointers to application procedures */
  bsf_BeginRead_fptr BeginReader;
  bsf_EndRead_fptr   EndReader;
  bsf_BC_fptr        BezierCurveReader;
  bsf_BSC_fptr       BSplineCurveReader;
  bsf_BP_fptr        BezierPatchReader;
  bsf_BSP_fptr       BSplinePatchReader;
  bsf_BSM_fptr       BSMeshReader;
  bsf_BSH_fptr       BSplineHoleReader;
  bsf_CPMark_fptr    CPMarkReader;
  bsf_Camera_fptr    CameraReader;
  bsf_Colour_fptr    ColourReader;
} bsf_UserReaders;

/* set up readers data structure to an empty state */
void bsf_ClearReaders ( bsf_UserReaders *readers );

/* register application's reading procedures */
void bsf_BeginReadingFuncd ( bsf_UserReaders *readers,
                             bsf_BeginRead_fptr BeginReader );
void bsf_EndReadingFuncd ( bsf_UserReaders *readers,
                           bsf_EndRead_fptr EndReader );
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
void bsf_CPMarkReadFunc ( bsf_UserReaders *readers,
                          bsf_CPMark_fptr CPMarkReader );
void bsf_CameraReadFuncd ( bsf_UserReaders *readers,
                           bsf_Camera_fptr CameraReader );
void bsf_ColourReadFuncd ( bsf_UserReaders *readers,
                           bsf_Colour_fptr ColourReader );

/* read the entire data file */
boolean bsf_ReadBSFiled ( const char *filename, bsf_UserReaders *readers );

/* auxiliary reading procedures, for internal use */
boolean _bsf_ReadCPMark ( bsf_UserReaders *readers, int maxnpoints );
boolean _bsf_ReadColour ( bsf_UserReaders *readers );

/* ////////////////////////////////////////////////////////////////////////// */
/* auxiliary, low-level reading procedures */
void bsf_GetNextSymbol ( void );

/* read partial data - for internal use */
boolean bsf_ReadIntNumber ( int *number );
boolean bsf_ReadDoubleNumber ( double *number );
boolean bsf_ReadPointd ( int maxcpdimen, double *point, int *cpdimen );
boolean bsf_ReadIdent ( int *ident );
int bsf_ReadPointsd ( int maxcpdimen, int maxnpoints,
                      double *points, int *cpdimen );
int bsf_ReadPointsMK ( int maxnpoints, unsigned int *mk );

boolean bsf_ReadSpaceDim ( int maxdim, int *spdimen );
boolean bsf_ReadCurveDegree ( int maxdeg, int *degree );
boolean bsf_ReadPatchDegree ( int maxdeg, int *udeg, int *vdeg );

boolean bsf_ReadKnotSequenced ( int maxlastknot, int *lastknot, double *knots,
                                boolean *closed );

boolean bsf_ReadCPMark ( int maxcp, unsigned int *mk );
boolean bsf_ReadCamera ( CameraRecd *Camera, int *ident );
boolean bsf_ReadColour ( point3d *colour );

/* the procedures with headers below are fit to read data if the application */
/* is sure about the file contents. The parameter readers should be NULL,    */
/* when called by an application. Otherwise use the high-level procedures,   */
/* whose prototypes were defined in the previous section of this file.       */
boolean bsf_OpenInputFile ( const char *filename );
void bsf_CloseInputFile ( void );
void bsf_PrintErrorLocation ( void );

boolean bsf_ReadBezierCurve4d ( int maxdeg, int *deg, point4d *cpoints,
                                int *spdimen, boolean *rational,
                                char *name, int *ident,
                                bsf_UserReaders *readers );
boolean bsf_ReadBSplineCurve4d ( int maxdeg, int maxlastknot, int maxncpoints,
                                 int *deg, int *lastknot, double *knots,
                                 boolean *closed, point4d *cpoints,
                                 int *spdimen, boolean *rational,
                                 char *name, int *ident,
                                 bsf_UserReaders *readers );
boolean bsf_ReadBezierPatch4d ( int maxdeg, int *udeg, int *vdeg,
                                int *pitch, point4d *cpoints,
                                int *spdimen, boolean *rational,
                                char *name, int *ident,
                                bsf_UserReaders *readers );
boolean bsf_ReadBSplinePatch4d ( int maxdeg, int maxlastknot, int maxncpoints,
                                 int *udeg, int *lastknotu, double *knotsu,
                                 int *vdeg, int *lastknotv, double *knotsv,
                                 boolean *closed_u, boolean *closed_v,
                                 int *pitch, point4d *cpoints,
                                 int *spdimen, boolean *rational,
                                 char *name, int *ident,
                                 bsf_UserReaders *readers );
boolean bsf_ReadBSMesh4d ( int maxnv, int maxnhe, int maxnfac,
                           int *degree,
                           int *nv, BSMvertex *mv, int *mvhei, point4d *vc,
                           int *nhe, BSMhalfedge *mhe,
                           int *nfac, BSMfacet *mfac, int *mfhei,
                           int *spdimen, boolean *rational,
                           char *name, int *ident,
                           bsf_UserReaders *readers );
boolean bsf_ReadBSplineHoled ( int maxk, int *hole_k, double *knots,
                               point2d *domain_cp, point4d *hole_cp,
                               int *spdimen, boolean *rational,
                               char *name, int *ident,
                               bsf_UserReaders *readers );

/* ////////////////////////////////////////////////////////////////////////// */
/* writing procedures */

/* the type definition below makes it possible to write additional attributes */
/* of objects (curvess nd surfaces), like colour. Applications are not        */
/* required to write such attributes and they may ignore them safely -- which */
/* requires some flexibility of the parameter lists and no need of changing   */
/* them in future.                                                            */
typedef boolean (*bsf_WriteAttr_fptr) ( void *userData );


boolean bsf_OpenOutputFile ( char *filename, boolean append );
void bsf_CloseOutputFile ( void );

void bsf_WriteCurrentIndentation ( void );
void bsf_WriteComment ( char *comment );

/* writing partial data */
void bsf_WriteDoubleNumber ( double x );
void bsf_WriteAltDoubleNumber ( double x, int nzf );
void bsf_WritePointd ( int cpdimen, const double *point );
void bsf_WritePointsd ( int cpdimen, int cols, int rows, int pitch,
                        const double *points );
void bsf_WriteSpaceDim ( int spdimen, boolean rational );
void bsf_WriteCurveDegree ( int degree );
void bsf_WritePatchDegree ( int udeg, int vdeg );
void bsf_WriteKnotSequenced ( int lastknot, const double *knots, boolean closed );
void bsf_WriteIdent ( int ident );

/* writing the entire geometric objects - curves and surfaces */
/* at the moment the only additional attribute written by the */
/* procedure passed as WriteAttr may be the colour, use the   */
/* bsf_WriteColour procedure for this. If no attributes need  */
/* to be written, the parameters WriteAttr and userData       */
/* may be NULL.                                               */
boolean bsf_WriteBezierCurved ( int spdimen, int cpdimen, boolean rational,
                                int deg, const double *cpoints,
                                const char *name, int ident,
                                bsf_WriteAttr_fptr WriteAttr, void *userData );
boolean bsf_WriteBSplineCurved ( int spdimen, int cpdimen, boolean rational,
                                 int deg, int lastknot, const double *knots,
                                 boolean closed, const double *cpoints,
                                 const char *name, int ident,
                                 bsf_WriteAttr_fptr WriteAttr, void *userData );
boolean bsf_WriteBezierPatchd ( int spdimen, int cpdimen, boolean rational,
                                int udeg, int vdeg, int pitch, const double *cpoints,
                                const char *name, int ident,
                                bsf_WriteAttr_fptr WriteAttr, void *userData );
boolean bsf_WriteBSplinePatchd ( int spdimen, int cpdimen, boolean rational,
                                 int udeg, int lastknotu, const double *knotsu,
                                 int vdeg, int lastknotv, const double *knotsv,
                                 boolean closed_u, boolean closed_v,
                                 int pitch, const double *cpoints,
                                 const char *name, int ident,
                                 bsf_WriteAttr_fptr WriteAttr, void *userData );
boolean bsf_WriteBSMeshd ( int spdimen, int cpdimen, boolean rational, int degree,
                           int nv, const BSMvertex *mv, const int *mvhei,
                           const double *vc,
                           int nhe, const BSMhalfedge *mhe,
                           int nfac, const BSMfacet *mfac, const int *mfhei,
                           const char *name, int ident,
                           bsf_WriteAttr_fptr WriteAttr, void *userData );
boolean bsf_WriteBSplineHoled ( int spdimen, int cpdimen, boolean rational,
                                int hole_k, const double *knots,
                                const point2d *domain_cp, const double *hole_cp,
                                const char *name, int ident,
                                bsf_WriteAttr_fptr WriteAttr, void *userData );

/* writing other objects and attributes */
void bsf_WritePointsMK ( int npoints, const unsigned int *mk );
boolean bsf_WriteColour ( point3d *colour );
boolean bsf_WriteCamera ( CameraRecd *Camera, int ident );

#ifdef __cplusplus
}
#endif

#endif /*BSFILE_H*/

