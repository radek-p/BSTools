
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libraybez library of C procedures -               */
/* ray tracing of Bezier patches                                         */ 

#ifndef CONST_  /* a dirty trick to suppress many compiler warnings */
#define CONST_ const
#endif

#ifndef RAYBEZD_H
#define RAYBEZD_H

#ifndef PKGEOM_H
#include "pkgeom.h"
#endif

#ifndef MULTIBS_H
#include "multibs.h"
#endif

#ifdef __cplusplus   
extern "C" {
#endif

typedef struct _BezPatchTreeVertexd {
  struct _BezPatchTreeVertexd *left, *right,  /* pointers to subtrees */
            *up;                              /* and up */
  point3d   *ctlpoints;                       /* pointer to array */
                                              /* of control points */
  vector3d  *normalvect;                      /* optional patch to describe */
                                              /* the normal vector */
  double    u0, u1, v0, v1;                   /* patch piece domain */
  Box3d     bbox;                             /* bounding box */
  point3d   pcent;                            /* patch central point */
  vector3d  nvcent;
  double    maxder;                           /* maximal derivative */
  short int level;
  char      divdir;                           /* 0 - divide u */
                                              /* 1 - divide v */
  boolean   leaf;
} BezPatchTreeVertexd, *BezPatchTreeVertexdp;

typedef struct {
  int                  object_id;  /* to be used by applications */
  unsigned char        n, m;       /* degree of the patch */
  unsigned int         cpsize;     /* size of array of control points */
  unsigned char        vn, vm;     /* degree of the normal vector patch */
  unsigned int         nvsize;     /* size of array of yhe normal vector patch */
                                   /* control points */
  BezPatchTreeVertexdp root;       /* root vertex representing the whole patch */
} BezPatchTreed, *BezPatchTreedp;


typedef struct _RBezPatchTreeVertexd {
  struct _RBezPatchTreeVertexd *left, *right,  /* pointers to subtrees */
            *up;                              /* and up */
  point4d   *ctlpoints;                       /* pointer to array */
                                              /* of control points */
  vector3d  *normalvect;                      /* optional patch to describe */
                                              /* the normal vector */
  double    u0, u1, v0, v1;                   /* patch piece domain */
  Box3d     bbox;                             /* bounding box */
  point3d   pcent;                            /* patch central point */
  vector3d  nvcent;
  double    maxder;                           /* maximal derivative */
  short int level;
  char      divdir;                           /* 0 - divide u */
                                              /* 1 - divide v */
  boolean   leaf;
} RBezPatchTreeVertexd, *RBezPatchTreeVertexdp;

typedef struct {
  int                  object_id;  /* to be used by applications */
  unsigned char        n, m;       /* degree of the patch */
  unsigned int         cpsize;     /* size of array of control points */
  unsigned char        vn, vm;     /* degree of the normal vector patch */
  unsigned int         nvsize;     /* size of array of yhe normal vector patch */
                                   /* control points */
  RBezPatchTreeVertexdp root;      /* root vertex representing the whole patch */
} RBezPatchTreed, *RBezPatchTreedp;


typedef struct _BezCurveTreeVertexd {
    struct _BezCurveTreeVertexd *left, *right,
                                *up;    /* pointers to subtrees and up */
    point3d *ctlpoints;                 /* pointer to array of control points */
    double  t0, t1;                     /* parameter range */
    Box3d   bbox;                       /* bounding box */
    point3d ccent;                      /* curve central point */
    double  maxder;
    short   level;                      /* subdivision level */
    boolean leaf;
    char    pad;
  } BezCurveTreeVertexd, *BezCurveTreeVertexdp;

typedef struct {
    int    object_id;           /* for applications */
    short  degree;              /* curve degree */
    short  cpsize;              /* size of array of control points */
    double ext;                 /* box extension */
    BezCurveTreeVertexdp root;  /* root vertex */
  } BezCurveTreed, *BezCurveTreedp;


typedef struct _RBezCurveTreeVertexd {
    struct _RBezCurveTreeVertexd *left, *right,
                                 *up;   /* pointers to subtrees and up */
    point4d *ctlpoints;                 /* pointer to array of control points */
    double  t0, t1;                     /* parameter range */
    Box3d   bbox;                       /* bounding box */
    point3d ccent;                      /* curve central point */
    double  maxder;
    short   level;                      /* subdivision level */
    boolean leaf;
    char    pad;
  } RBezCurveTreeVertexd, *RBezCurveTreeVertexdp;

typedef struct {
    int    object_id;            /* for applications */
    short  degree;               /* curve degree */
    short  cpsize;               /* size of array of control points */
    double ext;                  /* box extension */
    RBezCurveTreeVertexdp root;  /* root vertex */
  } RBezCurveTreed, *RBezCurveTreedp;


typedef struct {
  int      object_id;  /* to be used by applications */
  point3d  p;
  vector3d nv;
  double   u, v, t;
} RayObjectIntersd, *RayObjectIntersdp;


BezPatchTreedp
  rbez_NewBezPatchTreed ( int object_id,
                          unsigned char n, unsigned char m,
                          double u0, double u1, double v0, double v1,
                          CONST_ point3d *ctlpoints );

void rbez_DestroyBezPatchTreed ( BezPatchTreedp tree );

BezPatchTreeVertexdp
  rbez_GetBezLeftVertexd ( BezPatchTreedp tree,
                           BezPatchTreeVertexdp vertex );
BezPatchTreeVertexdp
  rbez_GetBezRightVertexd ( BezPatchTreedp tree,
                            BezPatchTreeVertexdp vertex );

int rbez_FindRayBezPatchIntersd ( BezPatchTreed *tree, ray3d *ray,
                                  int maxlevel, int maxinters,
                                  int *ninters, RayObjectIntersd *inters );


RBezPatchTreedp
  rbez_NewRBezPatchTreed ( int object_id,
                           unsigned char n, unsigned char m,
                           double u0, double u1, double v0, double v1,
                           CONST_ point4d *ctlpoints );

void rbez_DestroyRBezPatchTreed ( RBezPatchTreedp tree );

RBezPatchTreeVertexdp
  rbez_GetRBezLeftVertexd ( RBezPatchTreedp tree,
                            RBezPatchTreeVertexdp vertex );
RBezPatchTreeVertexdp
  rbez_GetRBezRightVertexd ( RBezPatchTreedp tree,
                             RBezPatchTreeVertexdp vertex );

int rbez_FindRayRBezPatchIntersd ( RBezPatchTreed *tree, ray3d *ray,
                                   int maxlevel, int maxinters,
                                   int *ninters, RayObjectIntersd *inters );


BezCurveTreedp rbez_NewBezCurveTreed ( int object_id, short degree,
                                       double t0, double t1, double ext,
                                       CONST_ point3d *ctlpoints );

void rbez_DestroyBezCurveTreed ( BezCurveTreedp tree );

BezCurveTreeVertexdp rbez_GetBezCurveLeftVertexd ( BezCurveTreedp tree,
                                                   BezCurveTreeVertexdp vertex );
BezCurveTreeVertexdp rbez_GetBezCurveRightVertexd ( BezCurveTreedp tree,
                                                    BezCurveTreeVertexdp vertex );

int rbez_FindRayBezcOffsetIntersd ( BezCurveTreedp tree, ray3d *ray,
                                    int maxlevel, int maxinters,
                                    int *ninters, RayObjectIntersd *inters );

RBezCurveTreedp rbez_NewRBezCurveTreed ( int object_id, short degree,
                                        double t0, double t1, double ext,
                                        CONST_ point4d *ctlpoints );

void rbez_DestroyRBezCurveTreed ( RBezCurveTreedp tree );

RBezCurveTreeVertexdp rbez_GetRBezCurveLeftVertexd ( RBezCurveTreedp tree,
                                                     RBezCurveTreeVertexdp vertex );
RBezCurveTreeVertexdp rbez_GetRBezCurveRightVertexd ( RBezCurveTreedp tree,
                                                      RBezCurveTreeVertexdp vertex );

int rbez_FindRayRBezcOffsetIntersd ( RBezCurveTreedp tree, ray3d *ray,
                                     int maxlevel, int maxinters,
                                     int *ninters, RayObjectIntersd *inters );

/* ////////////////////////////////////////////////////////////////////////// */
char rbez_TestRayBBoxd ( ray3d *ray, Box3d *box );

/* ////////////////////////////////////////////////////////////////////////// */
typedef struct {
    double u0, u1, v0, v1;
    double s0, s1, t0, t1;
    int    npoints;
  } rbiIntersArcd;

typedef void rbiArcOutd ( void *usrptr, rbiIntersArcd *arc, vector4d *ipt );

boolean rbi_FindRBezIntersectiond ( int n1, int m1, point4d *p1,
                                    int n2, int m2, point4d *p2,
                                    double epsilon, byte maxlevel,
                                    rbiArcOutd *outproc, void *usrptr );

/* ////////////////////////////////////////////////////////////////////////// */
boolean rbez_HomotopicClosedBSC3d ( int degree, int lastknot, double *knots,
                                    point3d *cpoints0, point3d *cpoints1,
                                    double *tfh, boolean *error );

/* ////////////////////////////////////////////////////////////////////////// */
#ifndef RAYBEZ_H
boolean raybez_InitMutex ( void );
void raybez_DestroyMutex ( void );
#endif

#ifdef __cplusplus
}
#endif

#endif

