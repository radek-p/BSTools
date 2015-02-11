
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
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

#define RAYBEZ_WHITE  0
#define RAYBEZ_GREY   1
#define RAYBEZ_BLACK  2

typedef struct _BezPatchTreeVertexd {
  struct _BezPatchTreeVertexd *left, *right,  /* pointers to subtrees */
                              *up;            /* and up */
  point3d                     *ctlpoints;     /* pointer to array */
                                              /* of control points */
  vector3d                    *nvcpoints;     /* optional patch to describe */
                                              /* the normal vector */
  double                      u0, u1, v0, v1; /* patch piece domain */
  Box3d                       bbox;           /* bounding box */
  point3d                     pcent;          /* patch central point */
  vector3d                    nvcent;         /* normal vector */
  double                      maxder;         /* maximal derivative length */
  unsigned char               level;          /* subdivision level */
  unsigned char               divdir;         /* 0 - divide u */
                                              /* 1 - divide v */
  unsigned char               vertex_colour;  /* for trimmed patches */
  unsigned char               tag;            /* for pthreads */
} BezPatchTreeVertexd, *BezPatchTreeVertexdp;

typedef struct {
  int                  object_id;  /* to be used by applications */
  unsigned char        n, m;       /* degree of the patch */
  unsigned char        nvn, nvm;   /* degree of the normal vector patch */
  unsigned int         cpsize;     /* size of array of control points */
  unsigned int         nvsize;     /* size of array of the normal vector patch */
                                   /* control points */
  BezPatchTreeVertexdp root;       /* root vertex representing the whole patch */
} BezPatchTreed, *BezPatchTreedp;


typedef struct _RBezPatchTreeVertexd {
  struct _RBezPatchTreeVertexd *left, *right,  /* pointers to subtrees */
                               *up;            /* and up */
  point4d                      *ctlpoints;     /* pointer to array */
                                               /* of control points */
  vector3d                     *nvcpoints;     /* optional patch to describe */
                                               /* the normal vector */
  double                       u0, u1, v0, v1; /* patch piece domain */
  Box3d                        bbox;           /* bounding box */
  point3d                      pcent;          /* patch central point */
  vector3d                     nvcent;         /* normal vector */
  double                       maxder;         /* maximal derivative length */
  unsigned char                level;          /* subdivision level */
  unsigned char                divdir;         /* 0 - divide u */
                                               /* 1 - divide v */
  unsigned char                vertex_colour;  /* for trimmed patches */
  unsigned char                tag;            /* for pthreads */
} RBezPatchTreeVertexd, *RBezPatchTreeVertexdp;

typedef struct {
  int                  object_id;  /* to be used by applications */
  unsigned char        n, m;       /* degree of the patch */
  unsigned char        vn, vm;     /* degree of the normal vector patch */
  unsigned int         cpsize;     /* size of array of control points */
  unsigned int         nvsize;     /* size of array of yhe normal vector patch */
                                   /* control points */
  RBezPatchTreeVertexdp root;      /* root vertex representing the whole patch */
} RBezPatchTreed, *RBezPatchTreedp;


typedef struct _BezCurveTreeVertexd {
    struct _BezCurveTreeVertexd *left, *right,
                                *up;        /* pointers to subtrees and up */
    point3d                     *ctlpoints; /* pointer to array of control points */
    double                      t0, t1;     /* parameter range */
    Box3d                       bbox;       /* bounding box */
    point3d                     ccent;      /* curve central point */
    double                      maxder;
    unsigned char               level;      /* subdivision level */
    unsigned char               tag;        /* for pthreads */
    unsigned short              pad;
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
                                 *up;        /* pointers to subtrees and up */
    point4d                      *ctlpoints; /* pointer to array of control points */
    double                       t0, t1;     /* parameter range */
    Box3d                        bbox;       /* bounding box */
    point3d                      ccent;      /* curve central point */
    double                       maxder;
    unsigned char                level;      /* subdivision level */
    unsigned char                tag;        /* for pthreads */
    unsigned short               pad;
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
  void     *extra_info;
} RayObjectIntersd, *RayObjectIntersdp;


/* ////////////////////////////////////////////////////////////////////////// */
void rbez_InitBBox3d ( Box3d *bbox, point3d *p );
void rbez_ExtendBBox3d ( Box3d *bbox, point3d *p );
void rbez_Extend2BBox3d ( Box3d *bbox, point3d *p1, point3d *p2 );

void rbez_FindCPBoundingBox3d ( int nrows, int ncols, int pitch,
                                point3d *cp, double eps, Box3d *bbox );
void rbez_FindCPBoundingBox3Rd ( int nrows, int ncols, int pitch,
                                 point4d *cp, double eps, Box3d *bbox );

void rbez_FindSumBBoxd ( Box3d *box1, Box3d *box2, Box3d *box );
boolean rbez_NarrowBBoxSumd ( Box3d *box1, Box3d *box2, Box3d *box );
boolean rbez_TestRayBBoxd ( ray3d *ray, Box3d *box );

/* ////////////////////////////////////////////////////////////////////////// */
/* the tree created by the two procedures below is built in the same way      */
BezPatchTreedp
  rbez_NewBezPatchTreed ( int object_id,
                          unsigned char n, unsigned char m,
                          double u0, double u1, double v0, double v1,
                          CONST_ point3d *ctlpoints );
BezPatchTreedp
    rbez_NewBSPatchTreed ( int object_id,
                   unsigned char n, unsigned int lknu, CONST_ double *knotsu,
                   unsigned char v, unsigned int lknv, CONST_ double *knotsv,
                   unsigned int pitch, CONST_ point3d *ctlpoints );

void rbez_DestroyBezPatchTreed ( BezPatchTreedp tree );

BezPatchTreeVertexdp
  rbez_GetBezLeftVertexd ( BezPatchTreedp tree,
                           BezPatchTreeVertexdp vertex );
BezPatchTreeVertexdp
  rbez_GetBezRightVertexd ( BezPatchTreedp tree,
                            BezPatchTreeVertexdp vertex );

int rbez_RayBezpWspSized ( int n, int m, int maxlevel );
int _rbez_FindRayBezPatchIntersd ( BezPatchTreed *tree, ray3d *ray,
                                   int maxlevel, int maxinters,
                                   int *ninters, RayObjectIntersd *inters,
                                   void *workspace );
int rbez_FindRayBezPatchIntersd ( BezPatchTreed *tree, ray3d *ray,
                                  int maxlevel, int maxinters,
                                  int *ninters, RayObjectIntersd *inters );

/* ////////////////////////////////////////////////////////////////////////// */
/* the tree created by the two procedures below is built in the same way      */
RBezPatchTreedp
  rbez_NewRBezPatchTreed ( int object_id,
                           unsigned char n, unsigned char m,
                           double u0, double u1, double v0, double v1,
                           CONST_ point4d *ctlpoints );
RBezPatchTreedp
    rbez_NewRBSPatchTreed ( int object_id,
                   unsigned char n, unsigned int lknu, CONST_ double *knotsu,
                   unsigned char v, unsigned int lknv, CONST_ double *knotsv,
                   unsigned int pitch, CONST_ point4d *ctlpoints );

void rbez_DestroyRBezPatchTreed ( RBezPatchTreedp tree );

RBezPatchTreeVertexdp
  rbez_GetRBezLeftVertexd ( RBezPatchTreedp tree,
                            RBezPatchTreeVertexdp vertex );
RBezPatchTreeVertexdp
  rbez_GetRBezRightVertexd ( RBezPatchTreedp tree,
                             RBezPatchTreeVertexdp vertex );

int rbez_RayRBezpWspSized ( int n, int m, int maxlevel );
int _rbez_FindRayRBezPatchIntersd ( RBezPatchTreed *tree, ray3d *ray,
                                    int maxlevel, int maxinters,
                                    int *ninters, RayObjectIntersd *inters,
                                    void *workspace );
int rbez_FindRayRBezPatchIntersd ( RBezPatchTreed *tree, ray3d *ray,
                                   int maxlevel, int maxinters,
                                   int *ninters, RayObjectIntersd *inters );

/* ////////////////////////////////////////////////////////////////////////// */
BezCurveTreedp rbez_NewBezCurveTreed ( int object_id, short degree,
                                       double t0, double t1,
                                       CONST_ point3d *ctlpoints, double ext );
BezCurveTreedp rbez_NewBSCurveTreed ( int object_id,
                                      short degree, int lastknot, CONST_ double *knots,
                                      CONST_ point3d *ctlpoints, double ext );

void rbez_DestroyBezCurveTreed ( BezCurveTreedp tree );

BezCurveTreeVertexdp rbez_GetBezCurveLeftVertexd ( BezCurveTreedp tree,
                                                   BezCurveTreeVertexdp vertex );
BezCurveTreeVertexdp rbez_GetBezCurveRightVertexd ( BezCurveTreedp tree,
                                                    BezCurveTreeVertexdp vertex );

int rbez_RayBezcOffsetWspSized ( int degree, int maxlevel );
int _rbez_FindRayBezcOffsetIntersd ( BezCurveTreedp tree, ray3d *ray,
                                     int maxlevel, int maxinters,
                                     int *ninters, RayObjectIntersd *inters,
                                     void *workspace );
int rbez_FindRayBezcOffsetIntersd ( BezCurveTreedp tree, ray3d *ray,
                                    int maxlevel, int maxinters,
                                    int *ninters, RayObjectIntersd *inters );

/* ////////////////////////////////////////////////////////////////////////// */
RBezCurveTreedp rbez_NewRBezCurveTreed ( int object_id, short degree,
                                        double t0, double t1,
                                        CONST_ point4d *ctlpoints, double ext );
RBezCurveTreedp rbez_NewRBSCurveTreed ( int object_id,
                                        short degree, int lastknot, CONST_ double *knots,
                                        CONST_ point4d *ctlpoints, double ext );

void rbez_DestroyRBezCurveTreed ( RBezCurveTreedp tree );

RBezCurveTreeVertexdp rbez_GetRBezCurveLeftVertexd ( RBezCurveTreedp tree,
                                                     RBezCurveTreeVertexdp vertex );
RBezCurveTreeVertexdp rbez_GetRBezCurveRightVertexd ( RBezCurveTreedp tree,
                                                      RBezCurveTreeVertexdp vertex );

int rbez_RayRBezcOffsetWspSized ( int degree, int maxlevel );
int _rbez_FindRayRBezcOffsetIntersd ( RBezCurveTreedp tree, ray3d *ray,
                                      int maxlevel, int maxinters,
                                      int *ninters, RayObjectIntersd *inters,
                                      void *workspace );
int rbez_FindRayRBezcOffsetIntersd ( RBezCurveTreedp tree, ray3d *ray,
                                     int maxlevel, int maxinters,
                                     int *ninters, RayObjectIntersd *inters );

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
boolean rbez_FindBezcHighlightPointsd ( int degree, point3d *cp,
                     double t0, double t1,
                     point3d *a, int maxlevel,
                     boolean (*out)(void *usrptr, double t, boolean singular),
                     void *usrptr );
boolean rbez_FindBezpHighlightPointsd ( int n, int m, point3d *cp,
                     double u0, double u1, double v0, double v1,
                     point3d *a, int maxlevel,
                     boolean (*out)(void *usrptr, point2d *uv, boolean singular),
                     void *usrptr );

boolean rbez_FindRBezcHighlightPointsd ( int degree, point4d *cp,
                     double t0, double t1,
                     point3d *a, int maxlevel,
                     boolean (*out)(void *usrptr, double t, boolean singular),
                     void *usrptr );
boolean rbez_FindRBezpHighlightPointsd ( int n, int m, point4d *cp,
                     double u0, double u1, double v0, double v1,
                     point3d *a, int maxlevel,
                     boolean (*out)(void *usrptr, point2d *uv, boolean singular),
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

