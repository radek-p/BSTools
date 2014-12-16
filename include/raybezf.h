
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libraybez library of C procedures -               */
/* ray tracing of Bezier patches                                         */ 

#ifndef CONST_  /* a dirty trick to suppress many compiler warnings */
#define CONST_ const
#endif

#ifndef RAYBEZF_H
#define RAYBEZF_H

#ifndef PKGEOM_H
#include "pkgeom.h"
#endif

#ifndef MULTIBS_H
#include "multibs.h"
#endif

#ifdef __cplusplus   
extern "C" {
#endif

typedef struct _BezPatchTreeVertexf {
  struct _BezPatchTreeVertexf *left, *right,  /* pointers to subtrees */
            *up;                              /* and up */
  point3f   *ctlpoints;                       /* pointer to array */
                                              /* of control points */
  vector3f  *normalvect;                      /* optional patch to describe */
                                              /* the normal vector */
  float     u0, u1, v0, v1;                   /* patch piece domain */
  Box3f     bbox;                             /* bounding box */
  point3f   pcent;                            /* patch central point */
  vector3f  nvcent;
  float     maxder;                           /* maximal derivative */
  short int level;
  char      divdir;                           /* 0 - divide u */
                                              /* 1 - divide v */
  boolean   leaf;
} BezPatchTreeVertexf, *BezPatchTreeVertexfp;

typedef struct {
  int                  object_id;  /* to be used by applications */
  unsigned char        n, m;       /* degree of the patch */
  unsigned int         cpsize;     /* size of array of control points */
  unsigned char        vn, vm;     /* degree of the normal vector patch */
  unsigned int         nvsize;     /* size of array of yhe normal vector patch */
                                   /* control points */
  BezPatchTreeVertexfp root;       /* root vertex representing the whole patch */
} BezPatchTreef, *BezPatchTreefp;


typedef struct _RBezPatchTreeVertexf {
  struct _RBezPatchTreeVertexf *left, *right,  /* pointers to subtrees */
            *up;                              /* and up */
  point4f   *ctlpoints;                       /* pointer to array */
                                              /* of control points */
  vector3f  *normalvect;                      /* optional patch to describe */
                                              /* the normal vector */
  float     u0, u1, v0, v1;                   /* patch piece domain */
  Box3f     bbox;                             /* bounding box */
  point3f   pcent;                            /* patch central point */
  vector3f  nvcent;
  float     maxder;                           /* maximal derivative */
  short int level;
  char      divdir;                           /* 0 - divide u */
                                              /* 1 - divide v */
  boolean   leaf;
} RBezPatchTreeVertexf, *RBezPatchTreeVertexfp;

typedef struct {
  int                  object_id;  /* to be used by applications */
  unsigned char        n, m;       /* degree of the patch */
  unsigned int         cpsize;     /* size of array of control points */
  unsigned char        vn, vm;     /* degree of the normal vector patch */
  unsigned int         nvsize;     /* size of array of yhe normal vector patch */
                                   /* control points */
  RBezPatchTreeVertexfp root;      /* root vertex representing the whole patch */
} RBezPatchTreef, *RBezPatchTreefp;


typedef struct _BezCurveTreeVertexf {
    struct _BezCurveTreeVertexf *left, *right,
                                *up;    /* pointers to subtrees and up */
    point3f *ctlpoints;                 /* pointer to array of control points */
    float   t0, t1;                     /* parameter range */
    Box3f   bbox;                       /* bounding box */
    point3f ccent;                      /* curve central point */
    float   maxder;
    short   level;                      /* subdivision level */
    boolean leaf;
    char    pad;
  } BezCurveTreeVertexf, *BezCurveTreeVertexfp;

typedef struct {
    int    object_id;           /* for applications */
    short  degree;              /* curve degree */
    short  cpsize;              /* size of array of control points */
    float ext;                  /* box extension */
    BezCurveTreeVertexfp root;  /* root vertex */
  } BezCurveTreef, *BezCurveTreefp;


typedef struct _RBezCurveTreeVertexf {
    struct _RBezCurveTreeVertexf *left, *right,
                                 *up;   /* pointers to subtrees and up */
    point4f *ctlpoints;                 /* pointer to array of control points */
    float   t0, t1;                     /* parameter range */
    Box3f   bbox;                       /* bounding box */
    point3f ccent;                      /* curve central point */
    float   maxder;
    short   level;                      /* subdivision level */
    boolean leaf;
    char    pad;
  } RBezCurveTreeVertexf, *RBezCurveTreeVertexfp;

typedef struct {
    int    object_id;            /* for applications */
    short  degree;               /* curve degree */
    short  cpsize;               /* size of array of control points */
    float  ext;                  /* box extension */
    RBezCurveTreeVertexfp root;  /* root vertex */
  } RBezCurveTreef, *RBezCurveTreefp;


typedef struct {
  int      object_id;  /* to be used by applications */
  point3f  p;
  vector3f nv;
  float    u, v, t;
} RayObjectIntersf, *RayObjectIntersfp;


/* ////////////////////////////////////////////////////////////////////////// */
void rbez_FindCPBoundingBox3f ( int nrows, int ncols, int pitch,
                                point3f *cp, float eps, Box3f *bbox );
void rbez_FindCPBoundingBox3Rf ( int nrows, int ncols, int pitch,
                                 point4f *cp, float eps, Box3f *bbox );

void rbez_FindSumBBoxf ( Box3f *box1, Box3f *box2, Box3f *box );
boolean rbez_NarrowBBoxSumf ( Box3f *box1, Box3f *box2, Box3f *box );
boolean rbez_TestRayBBoxf ( ray3f *ray, Box3f *box );

/* ////////////////////////////////////////////////////////////////////////// */
BezPatchTreefp
  rbez_NewBezPatchTreef ( int object_id,
                          unsigned char n, unsigned char m,
                          float u0, float u1, float v0, float v1,
                          CONST_ point3f *ctlpoints );

void rbez_DestroyBezPatchTreef ( BezPatchTreefp tree );

BezPatchTreeVertexfp
  rbez_GetBezLeftVertexf ( BezPatchTreefp tree,
                           BezPatchTreeVertexfp vertex );
BezPatchTreeVertexfp
  rbez_GetBezRightVertexf ( BezPatchTreefp tree,
                            BezPatchTreeVertexfp vertex );

int rbez_FindRayBezPatchIntersf ( BezPatchTreef *tree, ray3f *ray,
                                  int maxlevel, int maxinters,
                                  int *ninters, RayObjectIntersf *inters );


RBezPatchTreefp
  rbez_NewRBezPatchTreef ( int object_id,
                           unsigned char n, unsigned char m,
                           float u0, float u1, float v0, float v1,
                           CONST_ point4f *ctlpoints );

void rbez_DestroyRBezPatchTreef ( RBezPatchTreefp tree );

RBezPatchTreeVertexfp
  rbez_GetRBezLeftVertexf ( RBezPatchTreefp tree,
                            RBezPatchTreeVertexfp vertex );
RBezPatchTreeVertexfp
  rbez_GetRBezRightVertexf ( RBezPatchTreefp tree,
                             RBezPatchTreeVertexfp vertex );

int rbez_FindRayRBezPatchIntersf ( RBezPatchTreef *tree, ray3f *ray,
                                   int maxlevel, int maxinters,
                                   int *ninters, RayObjectIntersf *inters );


BezCurveTreefp rbez_NewBezCurveTreef ( int object_id, short degree,
                                       float t0, float t1, float ext,
                                       CONST_ point3f *ctlpoints );

void rbez_DestroyBezCurveTreef ( BezCurveTreefp tree );

BezCurveTreeVertexfp rbez_GetBezCurveLeftVertexf ( BezCurveTreefp tree,
                                                   BezCurveTreeVertexfp vertex );
BezCurveTreeVertexfp rbez_GetBezCurveRightVertexf ( BezCurveTreefp tree,
                                                    BezCurveTreeVertexfp vertex );

int rbez_FindRayBezcOffsetIntersf ( BezCurveTreefp tree, ray3f *ray,
                                    int maxlevel, int maxinters,
                                    int *ninters, RayObjectIntersf *inters );

RBezCurveTreefp rbez_NewRBezCurveTreef ( int object_id, short degree,
                                        float t0, float t1, float ext,
                                        CONST_ point4f *ctlpoints );

void rbez_DestroyRBezCurveTreef ( RBezCurveTreefp tree );

RBezCurveTreeVertexfp rbez_GetRBezCurveLeftVertexf ( RBezCurveTreefp tree,
                                                     RBezCurveTreeVertexfp vertex );
RBezCurveTreeVertexfp rbez_GetRBezCurveRightVertexf ( RBezCurveTreefp tree,
                                                      RBezCurveTreeVertexfp vertex );

int rbez_FindRayRBezcOffsetIntersf ( RBezCurveTreefp tree, ray3f *ray,
                                     int maxlevel, int maxinters,
                                     int *ninters, RayObjectIntersf *inters );

/* ////////////////////////////////////////////////////////////////////////// */
typedef struct {
    float u0, u1, v0, v1;
    float s0, s1, t0, t1;
    int   npoints;
  } rbiIntersArcf;

typedef void rbiArcOutf ( void *usrptr, rbiIntersArcf *arc, vector4f *ipt );

boolean rbi_FindRBezIntersectionf ( int n1, int m1, point4f *p1,
                                    int n2, int m2, point4f *p2,
                                    float epsilon, byte maxlevel,
                                    rbiArcOutf *outproc, void *usrptr );

/* ////////////////////////////////////////////////////////////////////////// */
boolean rbez_HomotopicClosedBSC3f ( int degree, int lastknot, float *knots,
                                    point3f *cpoints0, point3f *cpoints1,   
                                    float *tfh, boolean *error );

/* ////////////////////////////////////////////////////////////////////////// */
boolean rbez_FindBezcHighlightPointsf ( int degree, point3f *cp,
                     float t0, float t1,
                     point3f *a, int maxlevel,
                     boolean (*out)(void *usrptr, float t, boolean singular),
                     void *usrptr );
boolean rbez_FindBezpHighlightPointsf ( int n, int m, point3f *cp,
                     float u0, float u1, float v0, float v1,
                     point3f *a, int maxlevel,
                     boolean (*out)(void *usrptr, point2f *uv, boolean singular),
                     void *usrptr );

boolean rbez_FindRBezcHighlightPointsf ( int degree, point4f *cp,
                     float t0, float t1,
                     point3f *a, int maxlevel,
                     boolean (*out)(void *usrptr, float t, boolean singular),
                     void *usrptr );
boolean rbez_FindRBezpHighlightPointsf ( int n, int m, point4f *cp,
                     float u0, float u1, float v0, float v1,   
                     point3f *a, int maxlevel,
                     boolean (*out)(void *usrptr, point2f *uv, boolean singular),
                     void *usrptr );

/* ////////////////////////////////////////////////////////////////////////// */
#ifndef RAYBEZ_H
boolean raybez_InitMutex ( void );
void raybez_DestroyMutex ( void );
#endif

#ifdef __cplusplus
}
#endif

#endif

