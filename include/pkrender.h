
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */  
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef PKRENDER_H
#define PKRENDER_H

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif
#ifndef PKVTHREADS_H
#include "pkvthreads.h"
#endif
#ifndef PKGEOM_H
#include "pkgeom.h"
#endif
#ifndef PKNUM_H
#include "pknum.h"
#endif
#ifndef CAMERA_H
#include "camera.h"
#endif
#ifndef MULTIBS_H
#include "multibs.h"
#endif
#ifndef RAYBEZ_H
#include "raybez.h"
#endif

/* sizes of stacks and scratch memory pools for threads */
#define PKRENDER_STACK       4*1048576
#define PKRENDER_SCRATCHMEM 16*1048576

#define MAXPLEVEL      15
#define MAXCLEVEL      24
#define MAXINTERS      32

#define R_NLIGHTS       4

/* scaling factor range for discontinuous palette */
#define R_MINDFSF        0.1
#define R_MAXDFSF     1000.0

#define MIN_T 1.0e-5

/* shape function type */
#define shapefunc_NONE         0
#define shapefunc_GAUSSIAN     1
#define shapefunc_MEAN         2
#define shapefunc_REFLECTION   3
#define shapefunc_HIGHLIGHT    4
#define shapefunc_LAMBERTIAN   5
#define shapefunc_CROSSECTIONS 6

/* object types, more to be added */
#define obj_TRIANGLE  1
#define obj_BSPATCH   2
#define obj_RBSPATCH  3
#define obj_BEZCURVE  4
#define obj_RBEZCURVE 5

typedef struct {
    Box3d           bbox;
    point3d         p0;
    vector3d        a1, a2, n;
  } triangle_data;

typedef struct {
    char            type;  /* == obj_TRIANGLE */
    double          colour[3];
    triangle_data   *trdata;
  } renderobj_triangle;

typedef struct {
    char            type;  /* == obj_BSPATCH */
    double          colour[3];
    BezPatchTreedp  ptree;
  } renderobj_bspatch;

typedef struct {
    char            type;  /* == obj_RBSPATCH */
    double          colour[3];
    RBezPatchTreedp ptree;
  } renderobj_rbspatch;

typedef struct {
    char            type;  /* == obj_BEZCURVE */
    double          colour[3];
    BezCurveTreedp  ctree;
  } renderobj_bezcurve;

typedef struct {
    char            type;  /* == obj_RBEZCURVE */
    double          colour[3];
    RBezCurveTreedp ctree;
  } renderobj_rbezcurve;

typedef union {
    char type;
    renderobj_triangle  triang;
    renderobj_bspatch   bsp;
    renderobj_rbspatch  rbsp;
    renderobj_bezcurve  bezc;
    renderobj_rbezcurve rbezc;
  } renderobj;

/* binary tree of object hierarchy */
typedef struct tnode {
    Box3d        bbox;
    struct tnode *left, *right;
    int          k;
  } rendertnode;

typedef struct {
    short       maxwidth, maxheight;
    void        (*SetPixel)( void *private_data,
                             short x, short y, byte r, byte g, byte b );
    CameraRecd  CPos;
    void        *private_data;
    byte        *imagedata;
      /* array of objects */
    int         obj_tab_length, nobjects;
    renderobj   *obj_tab;
    rendertnode *root;

    boolean     swShadows, swAntialias;
    vector3d    lightdir[R_NLIGHTS];
    double      lightint[R_NLIGHTS+1];

      /* switches selecting the shape function for visualisation */
    short int   c_shape_func, d_shape_func;
      /* shape function data */
    vector3d    sectiondir;
    point3d     rp0;
    vector3d    rv1, rv2;
    point3d     hp0;
    vector3d    hv1, hv2;
    double      minsf, maxsf;
      /* scaling factor for rendering the shape function */
    double      dfsf;
    double      cfrange[2];
      /* number of rendering threads */
    short       nthr;
      /* rendering in progress */
    boolean     RenderingIsOn;
    short       y;
    byte        *aabuf, *aaline[7];
    clock_t     tic, ticks1, ticks2;
  } pkRenderer;


/* ///////////////////////////////////////////////////////////////////////// */
/* to be called at the beginning and end of the application */
boolean RendInit ( pkRenderer *rend, void *private_data, short nthr,
                   short maxwidth, short maxheight,
                   void (*SetPixel)( void *private_data,
                   short x, short y, byte r, byte g, byte b ) );
void RendDestroy ( pkRenderer *rend );

/* to be called before rendering a new scene */
boolean RendReset ( pkRenderer *rend );

boolean RendEnterTriangle3d ( pkRenderer *rend,
                              point3d *p0, point3d *p1, point3d *p2,
                              double *colour );
boolean RendEnterBezPatch3d ( pkRenderer *rend,
                              int n, int m, const point3d *cp, double *colour );
boolean RendEnterBSPatch3d ( pkRenderer *rend,
                             int n, int lknu, const double *knu,
                             int m, int lknv, const double *knv,
                             const point3d *cp, double *colour );
boolean RendEnterBezPatch3Rd ( pkRenderer *rend,
                               int n, int m, const point4d *cp, double *colour );
boolean RendEnterBSPatch3Rd ( pkRenderer *rend,
                              int n, int lknu, const double *knu,
                              int m, int lknv, const double *knv,
                              const point4d *cp, double *colour );
boolean RendEnterBezCurve3d ( pkRenderer *rend,
                              int n, const point3d *cp, double r, double *colour );
boolean RendEnterBSCurve3d ( pkRenderer *rend, int n, int lkn, const double *kn,
                             const point3d *cp, double r, double *colour );
boolean RendEnterBezCurve3Rd ( pkRenderer *rend,
                               int n, const point4d *cp, double r, double *colour );
boolean RendEnterBSCurve3Rd ( pkRenderer *rend, int n, int lkn, const double *kn,
                              const point4d *cp, double r, double *colour );

boolean RendEnterCamerad ( pkRenderer *rend, CameraRecd *CPos );
void RendEnterLightsd ( pkRenderer *rend,
                        int nlights, const vector3d *lightdir,
                        const double *lightint );
void RendEnterReflectionLinesFramed ( pkRenderer *rend, point3d rf[3] );
void RendEnterHighlightLinesFramed ( pkRenderer *rend, point3d hf[3] );
void RendEnterSectionPlanesNormald ( pkRenderer *rend, vector3d *spn );

/* actual rendering */
boolean RendBegin ( pkRenderer *rend, boolean antialias, short nthr );
boolean RendRestart ( pkRenderer *rend  );
void InitRenderingAA ( pkRenderer *rend  );
int RenderLineA ( pkRenderer *rend  );
int RenderLineAA ( pkRenderer *rend  );
int RenderLine ( pkRenderer *rend  );

#endif /*PKRENDER_H*/

