
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>   
#include <stdio.h>
#include <math.h>  
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/gl.h>  
#include <GL/glu.h> 
#include <GL/glx.h>

#include "pkvaria.h" 
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "bsmesh.h"
#include "g2blendingd.h"
#include "egholed.h"
#include "mengerc.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

#include "editor.h"
#include "editor_bsm.h"


void DrawAVertex ( int cpdimen, int spdimen, const double *cp )
{
  switch ( cpdimen ) {
case 2:
    glVertex2dv ( cp );
    break;
case 3:
    if ( spdimen == 2 )
      glVertex4d ( cp[0], cp[1], 0.0, cp[2] );
    else
      glVertex3dv ( cp );
    break;
case 4:
    glVertex4dv ( cp );
    break;
default:
    break;
  }
} /*DrawAVertex*/

void DrawAPolyline ( int cpdimen, int spdimen, int np, const double *cp )
{
  int i, k;

  glBegin ( GL_LINE_STRIP );
    switch ( cpdimen ) {
  case 2:
      for ( i = k = 0;  i < np;  i++, k += 2 )
        glVertex2dv ( &cp[k] );
      break;
  case 3:
      if ( spdimen == 2 ) {
        for ( i = k = 0;  i < np;  i++, k += 3 )
          glVertex4d ( cp[k], cp[k+1], 0.0, cp[k+2] );
      }
      else {
        for ( i = k = 0;  i < np;  i++, k += 3 )
          glVertex3dv ( &cp[k] );
      }
      break;
  case 4:
      for ( i = k = 0;  i < np;  i++, k += 4 )
        glVertex4dv ( &cp[k] );
      break;
  default:
      break;
    }
  glEnd ();
} /*DrawAPolyline*/

void DrawLineSegments ( int cpdimen, int spdimen, int np,
                        const double *cp, const double *wp )
{
  int i, k;

  glBegin ( GL_LINES );
    switch ( cpdimen ) {
  case 2:
      for ( i = k = 0;  i < np;  i++, k += 2 ) {
        glVertex2dv ( &cp[k] );
        glVertex2dv ( &wp[k] );
      }
      break;
  case 3:
      if ( spdimen == 2 ) {
        for ( i = k = 0;  i < np;  i++, k += 3 ) {
          glVertex4d ( cp[k], cp[k+1], 0.0, cp[k+2] );
          glVertex4d ( wp[k], wp[k+1], 0.0, wp[k+2] );
        }
      }
      else {
        for ( i = k = 0;  i < np;  i++, k += 3 ) {
          glVertex3dv ( &cp[k] );
          glVertex3dv ( &wp[k] );
        }
      }
      break;
  case 4:
      for ( i = k = 0;  i < np;  i++, k += 4 ) {
        glVertex4dv ( &cp[k] );
        glVertex4dv ( &wp[k] );
      }
      break;
  default:
      break;
    }
  glEnd ();
} /*DrawLineSegments*/

void DrawARectNet ( int cpdimen, int spdimen, int ncols, int nrows,
                    int pitch, const double *cp )
{
  int i, j, k;

  for ( i = k = 0;  i < ncols;  i++, k = i*pitch ) {
    glBegin ( GL_LINE_STRIP );
      switch ( cpdimen ) {
    case 2:
        for ( j = 0;  j < nrows;  j++, k += 2 )
          glVertex2dv ( &cp[k] );
        break;
    case 3:
        if ( spdimen == 2 ) {
          for ( j = 0;  j < nrows;  j++, k += 3 )
            glVertex4d ( cp[k], cp[k+1], 0.0, cp[k+2] );
        }
        else {
          for ( j = 0;  j < nrows;  j++, k += 3 )
            glVertex3dv ( &cp[k] );
        }
        break;
    case 4:
        for ( j = 0;  j < nrows;  j++, k += 4 )
          glVertex4dv ( &cp[k] );
        break;
    default:
        return;
      }
    glEnd ();
  }
  for ( j = 0; j < nrows; j++ ) {
    glBegin ( GL_LINE_STRIP );
      switch ( cpdimen ) {
    case 2:
        for ( i = 0, k = 2*j;  i < ncols;  i++, k += pitch )
          glVertex2dv ( &cp[k] );
        break;
    case 3:
        if ( spdimen == 2 ) {
          for ( i = 0, k = 3*j;  i < ncols;  i++, k += pitch )
            glVertex4d ( cp[k], cp[k+1], 0.0, cp[k+2] );
        }
        else {
          for ( i = 0, k = 3*j;  i < ncols;  i++, k += pitch )
            glVertex3dv ( &cp[k] );
        }
        break;
    case 4:
        for ( i = 0, k = 4*j;  i < ncols;  i++, k += pitch )
          glVertex4dv ( &cp[k] );
        break;
      }
    glEnd ();
  }
} /*DrawARectNet*/

void DrawCPoints ( int cpdimen, int spdimen, int np, const double *cp,
                   const byte *mkcp, byte mask1, byte mask2 )
{
  int  i, k;
  byte m;
#define IF_DRAWIT \
  m = mkcp[i];  if ( m & mask1 ) m |= MASK_MARKED; \
  if ( m == mask2 || (m & MASK_CP_SPECIAL) == mask2 )

  glPointSize ( 3.0 );
  glBegin ( GL_POINTS );
    switch ( cpdimen ) {
  case 2:
      for ( i = k = 0;  i < np;  i++, k += 2 )
        { IF_DRAWIT glVertex2dv ( &cp[k] ); }
      break;
  case 3:
      if ( spdimen == 2 ) {
        for ( i = k = 0;  i < np;  i++, k += 3 )
          { IF_DRAWIT glVertex4d ( cp[k], cp[k+1], 0.0, cp[k+2] ); }
      }
      else {
        for ( i = k = 0;  i < np;  i++, k += 3 )
          { IF_DRAWIT glVertex3dv ( &cp[k] ); }
      }
      break;
  case 4:
      for ( i = k = 0;  i < np;  i++, k += 4 )
        { IF_DRAWIT glVertex4dv ( &cp[k] ); }
      break;
  default:
      break;
    }
  glEnd ();
#undef IF_DRAWIT
} /*DrawCPoints*/

void DrawPoints ( int cpdimen, int spdimen, int np, const double *cp )
{
  int  i, k;

  glPointSize ( 3.0 );
  glBegin ( GL_POINTS );
    switch ( cpdimen ) {
  case 2:
      for ( i = k = 0;  i < np;  i++, k += 2 )
        glVertex2dv ( &cp[k] );
      break;
  case 3:
      if ( spdimen == 2 ) {
        for ( i = k = 0;  i < np;  i++, k += 3 )
          glVertex4d ( cp[k], cp[k+1], 0.0, cp[k+2] );
      }
      else {
        for ( i = k = 0;  i < np;  i++, k += 3 )
          glVertex3dv ( &cp[k] );
      }
      break;
  case 4:
      for ( i = k = 0;  i < np;  i++, k += 4 )
        glVertex4dv ( &cp[k] );
      break;
  default:
      break;
    }
  glEnd ();
} /*DrawPoints*/

void DrawBezierCurve ( int cpdimen, int spdimen, int deg, double *cp, int dd )
{
  int    i;
  double t, p[4];

  glBegin ( GL_LINE_STRIP );
  switch ( cpdimen ) {
case 2:
    for ( i = 0;  i <= dd;  i++ ) {
      t = (double)i/(double)dd;
      mbs_BCHornerC2d ( deg, cp, t, p );
      glVertex2dv ( p );
    }
    break;
case 3:
    if ( spdimen == 2 ) {
      for ( i = 0;  i <= dd;  i++ ) {
        t = (double)i/(double)dd;
        mbs_BCHornerC3d ( deg, cp, t, p );
        glVertex4d ( p[0], p[1], 0.0, p[2] );
      }
    }
    else {
      for ( i = 0;  i <= dd;  i++ ) {
        t = (double)i/(double)dd;
        mbs_BCHornerC3d ( deg, cp, t, p );
        glVertex3dv ( p );
      }
    }
    break;
case 4:
    for ( i = 0;  i <= dd;  i++ ) {
      t = (double)i/(double)dd;
      mbs_BCHornerC4d ( deg, cp, t, p );
      glVertex4dv ( p );
    }
    break;
default:
    break;
  }
  glEnd ();
} /*DrawBezierCurve*/

void DrawBSplineCurve ( int cpdimen, int spdimen, int deg, int lkn, double *knots,
                        double *cp, int dd )
{
  void   *sp;
  double *bcp;
  int    i, j, k, l;
  double t, p[4];

  sp = pkv_GetScratchMemTop ();
  l = mbs_NumKnotIntervalsd ( deg, lkn, knots );
  if ( l > 0 ) {
    bcp = pkv_GetScratchMemd ( l*(deg+1)*cpdimen );
    if ( bcp ) {
      mbs_multiBSCurvesToBezd ( cpdimen, 1, deg, lkn, knots, 0, cp,
                                &l, &k, NULL, 0, bcp );
      glBegin ( GL_LINE_STRIP );
      switch ( cpdimen ) {
    case 2:
        glVertex2dv ( bcp );
        for ( j = k = 0;  j < l;  j++, k += 2*(deg+1) )
          for ( i = 1;  i <= dd;  i++ ) {
            t = (double)i/(double)dd;
            mbs_BCHornerC2d ( deg, &bcp[k], t, p );
            glVertex2dv ( p );
          }
        break;
    case 3:
        if ( spdimen == 2 ) {
          glVertex4d ( bcp[0], bcp[1], 0.0, bcp[2] );
          for ( j = k = 0;  j < l;  j++, k += 3*(deg+1) )
            for ( i = 1;  i <= dd;  i++ ) {
              t = (double)i/(double)dd;
              mbs_BCHornerC3d ( deg, &bcp[k], t, p );
              glVertex4d ( p[0], p[1], 0.0, p[2] );
            }
        }
        else {
          glVertex3dv ( bcp );
          for ( j = k = 0;  j < l;  j++, k += 3*(deg+1) )
            for ( i = 1;  i <= dd;  i++ ) {
              t = (double)i/(double)dd;
              mbs_BCHornerC3d ( deg, &bcp[k], t, p );
              glVertex3dv ( p );
            }
        }
        break;
    case 4:
        glVertex4dv ( bcp );
        for ( j = k = 0;  j < l;  j++, k += 4*(deg+1) )
          for ( i = 1;  i <= dd;  i++ ) {
            t = (double)i/(double)dd;
            mbs_BCHornerC4d ( deg, &bcp[k], t, p );
            glVertex4dv ( p );
          }
        break;
    default:
        break;
      }
      glEnd ();
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*DrawBSplineCurve*/

void DrawBezierPatchWF ( int cpdimen, int spdimen, int degu, int degv,
                         int pitch, const double *cpoints,
                         int densu, int densv, int denscu, int denscv,
                         boolean firstu, boolean lastu,
                         boolean firstv, boolean lastv )
{
  void   *sp;
  double *cc, *dd;
  double t;
  int    i, j, size, fu, lu, fv, lv;

  sp = pkv_GetScratchMemTop ();
  size = (max ( degu, degv ) + 1)*cpdimen;
  cc = pkv_GetScratchMemd ( 2*size );
  if ( cc ) {
    dd = &cc[size];
    if ( firstu ) fu = 0;
             else fu = 1;
    if ( lastu ) lu = densu;
            else lu = densu-1;
    for ( i = fu; i <= lu; i++ ) {
      t = (double)i/(double)densu;
      if ( pitch == cpdimen*(degv+1) )
        mbs_multiBCHornerd ( degu, 1, cpdimen*(degv+1), 0, cpoints, t, cc );
      else {
              /* the columns are separated by gaps in the array */
        for ( j = 0; j <= degv; j++ ) {
          pkv_Selectd ( degu+1, cpdimen, pitch, cpdimen, &cpoints[j*cpdimen], dd );
          mbs_multiBCHornerd ( degu, 1, cpdimen, 0, dd, t, &cc[j*cpdimen] );
        }
      }
      DrawBezierCurve ( cpdimen, spdimen, degv, cc, denscv );
    }
    if ( firstv ) fv = 0;
             else fv = 1;
    if ( lastv ) lv = densv;
            else lv = densv-1;
    for ( i = fv; i <= lv; i++ ) {
      t = (double)i/(double)densv;
      mbs_multiBCHornerd ( degv, degu+1, cpdimen, pitch, cpoints, t, cc );
      DrawBezierCurve ( cpdimen, spdimen, degu, cc, denscu );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*DrawBezierPatchWF*/

void DrawBSplinePatchWF ( int cpdimen, int spdimen,
                          int degu, int lknu, const double *knu,
                          int degv, int lknv, const double *knv,
                          int pitch, const double *cpoints,
                          int densu, int densv, boolean firstu, boolean lastu,
                          boolean firstv, boolean lastv )
{
  void   *sp;
  int    ku, kv, bpitch;
  double *bcp;
  int    i, j;

  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsd ( degu, lknu, knu );
  kv = mbs_NumKnotIntervalsd ( degv, lknv, knv );
  bpitch = (degv+1)*kv*cpdimen;
  bcp = pkv_GetScratchMemd ( bpitch*ku*(degu+1) );
  if ( bcp ) {
    mbs_BSPatchToBezd ( cpdimen, degu, lknu, knu, degv, lknv, knv,
                        pitch, cpoints, &ku, NULL, NULL, &kv, NULL, NULL,
                        bpitch, bcp );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ )
        DrawBezierPatchWF ( cpdimen, spdimen, degu, degv, bpitch,
                            &bcp[i*(degu+1)*bpitch + j*(degv+1)*cpdimen],
                            densu, densv, 8*densu, 8*densv,
                            firstu || i > 0, lastu && i == ku-1,
                            firstv || j > 0, lastv && j == kv-1 );
  }
  pkv_SetScratchMemTop ( sp );
} /*DrawBSplinePatchWF*/

void GeomObjectDrawMeshEdges ( int cpdimen, int spdimen,
                               int nv, BSMvertex *mv, double *mvc, int *mvhei,
                               int nhe, BSMhalfedge *mhe,
                               int nfac, BSMfacet *mfac, int *mfhei,
                               byte edgemask )
{
  int     i, j;
  boolean ok, inner, boundary;

  glBegin ( GL_LINES );
    boundary = edgemask & MASK_HE_BOUNDARY;
    inner    = edgemask & MASK_HE_INNER;
    switch ( cpdimen ) {
  case 2:
      for ( i = 0; i < nhe; i++ ) {
        if ( mhe[i].otherhalf < 0 )
          ok = boundary;
        else if ( mhe[i].v0 < mhe[i].v1 )
          ok = inner;
        else
          ok = false;
        if ( ok ) {
          glVertex2dv ( &mvc[2*mhe[i].v0] );
          glVertex2dv ( &mvc[2*mhe[i].v1] );
        }
      }
      break;
  case 3:
      if ( spdimen == 2 ) {
        for ( i = 0; i < nhe; i++ ) {
          if ( mhe[i].otherhalf < 0 )
            ok = boundary;
          else if ( mhe[i].v0 < mhe[i].v1 )
            ok = inner;
          else
            ok = false;
          if ( ok ) {
            j = 3*mhe[i].v0;
            glVertex4d ( mvc[j], mvc[j+1], 0.0, mvc[j+2] );
            j = 3*mhe[i].v1;
            glVertex4d ( mvc[j], mvc[j+1], 0.0, mvc[j+2] );
          }
        }
      }
      else {
        for ( i = 0; i < nhe; i++ ) {
          if ( mhe[i].otherhalf < 0 )
            ok = boundary;
          else if ( mhe[i].v0 < mhe[i].v1 )
            ok = inner;
          else
            ok = false;
          if ( ok ) {
            glVertex3dv ( &mvc[3*mhe[i].v0] );
            glVertex3dv ( &mvc[3*mhe[i].v1] );
          }
        }
      }
      break;
  case 4:
      for ( i = 0; i < nhe; i++ ) {
        if ( mhe[i].otherhalf < 0 )
          ok = boundary;
        else if ( mhe[i].v0 < mhe[i].v1 )
          ok = inner;
        else
          ok = false;
        if ( ok ) {
          glVertex4dv ( &mvc[4*mhe[i].v0] );
          glVertex4dv ( &mvc[4*mhe[i].v1] );
        }
      }
      break;
  default:
      break;
    }
  glEnd ();
} /*GeomObjectDrawMeshEdges*/

void GeomObjectDrawMeshFacet ( int cpdimen, int spdimen,
                               int nv, BSMvertex *mv, double *mvc, int *mvhei,
                               int nhe, BSMhalfedge *mhe,
                               int nfac, BSMfacet *mfac, int *mfhei,
                               int facetnum )
{
#define SCF 0.75
  void   *sp;
  int    d, fhe, i, k;
  double *fv;

  if ( facetnum < 0 || facetnum > nfac )
    return;

  sp = pkv_GetScratchMemTop ();
  d = mfac[facetnum].degree;
  fv = pkv_GetScratchMemd ( (d+1)*cpdimen );
  if ( !fv )
    goto way_out;
        /* compute the facet central point */
  fhe = mfac[facetnum].firsthalfedge;
  memset ( fv, 0, cpdimen*sizeof(double) );
  for ( i = 0; i < d; i++ ) {
    memcpy ( &fv[(i+1)*cpdimen], &mvc[mhe[mfhei[fhe+i]].v0*cpdimen],
             cpdimen*sizeof(double) );
    pkn_AddMatrixd ( 1, cpdimen, 0, fv, 0, &fv[(i+1)*cpdimen], 0, fv );
  }
  pkn_MultMatrixNumd ( 1, cpdimen, 0, fv, 1.0/(double)d, 0, fv );
  pkn_MatrixLinCombd ( d, cpdimen, 0, fv, (1.0-SCF), cpdimen, &fv[cpdimen], SCF,
                       cpdimen, &fv[cpdimen] );
  glBegin ( GL_TRIANGLE_FAN );
    switch ( cpdimen ) {
  case 2:
      for ( i = k = 0;  i <= d;  i++, k += 2 )
        glVertex2dv ( &fv[k] );
      glVertex2dv ( &fv[cpdimen] );
      break;
  case 3:
      if ( spdimen == 2 ) {
        for ( i = k = 0;  i <= d;  i++, k += 3 )
          glVertex4d ( fv[k], fv[k+1], 0.0, fv[k+2] );
        glVertex4d ( fv[3], fv[4], 0.0, fv[5] );
      }
      else if ( spdimen == 3 ) {
        for ( i = k = 0;  i <= d;  i++, k += 3 )
          glVertex3dv ( &fv[k] );
        glVertex3dv ( &fv[cpdimen] );
      }
      break;
  case 4:
      for ( i = k = 0;  i <= d;  i++, k += 4 )
        glVertex4dv ( &fv[k] );
      glVertex4dv ( &fv[cpdimen] );
      break;
  default:
      break;
    }
  glEnd ();

way_out:
  pkv_SetScratchMemTop ( sp );
#undef SCF
} /*GeomObjectDrawMeshFacet*/

void GeomObjectDrawMeshHalfedge ( int cpdimen, int spdimen,
                                  int nv, BSMvertex *mv, double *mvc, int *mvhei,
                                  int nhe, BSMhalfedge *mhe,
                                  int nfac, BSMfacet *mfac, int *mfhei,
                                  int henum )
{
  int v0, v1;

  glLineWidth ( 3.0 );
  glBegin ( GL_LINES );
  v0 = mhe[henum].v0;
  v1 = mhe[henum].v1;
  switch ( cpdimen ) {
case 2:
    glVertex2dv ( &mvc[2*v0] );
    glVertex2dv ( &mvc[2*v1] );
    break;
case 3:
    if ( spdimen == 3 ) {
      glVertex3dv ( &mvc[3*v0] );
      glVertex3dv ( &mvc[3*v1] );
    }
    else {
      glVertex4d ( mvc[3*v0], mvc[3*v0+1], 0.0, mvc[3*v0+2] );
      glVertex4d ( mvc[3*v1], mvc[3*v1+1], 0.0, mvc[3*v1+2] );
    }
    break;
case 4:
    glVertex4dv ( &mvc[4*v0] );
    glVertex4dv ( &mvc[4*v1] );
    break;
default:
    break;
  }
  glEnd ();
  glLineWidth ( 1.0 );
} /*GeomObjectDrawMeshHalfedge*/

void GeomObjectDrawMeshHalfedges ( int cpdimen, int spdimen,
                                   int nv, BSMvertex *mv, double *mvc, int *mvhei,
                                   int nhe, BSMhalfedge *mhe,
                                   int nfac, BSMfacet *mfac, int *mfhei,
                                   byte edgemask,
                                   byte *hemark, byte hemask, boolean inv,
                                   boolean half0, boolean half1 )
{
  int     i, j, v0, v1;
  boolean ok, inner, boundary;
  double  midpoint[4];

#define CHECKIT1 \
  { \
    if ( inv ) { \
      if ( mhe[i].otherhalf < 0 ) \
        ok = boundary && !(hemark[i] & hemask); \
      else \
        ok = inner && !(hemark[i] & hemask); \
    } \
    else { \
      if ( mhe[i].otherhalf < 0 ) \
        ok = boundary && (hemark[i] & hemask); \
      else \
        ok = inner && (hemark[i] & hemask); \
    } \
  }

#define CHECKIT2 \
  { \
    if ( inv ) { \
      if ( mhe[i].otherhalf < 0 ) \
        ok = boundary && !(hemark[i] & hemask); \
      else if ( mhe[i].v0 < mhe[i].v1 ) \
        ok = inner && !((hemark[i] | hemark[mhe[i].otherhalf]) & hemask); \
      else \
        ok = false; \
    } \
    else { \
      if ( mhe[i].otherhalf < 0 ) \
        ok = boundary && (hemark[i] & hemask); \
      else if ( mhe[i].v0 < mhe[i].v1 ) \
        ok = inner && ((hemark[i] | hemark[mhe[i].otherhalf]) & hemask); \
      else \
        ok = false; \
    } \
  }

  glBegin ( GL_LINES );
    boundary = edgemask & MASK_HE_BOUNDARY;
    inner    = edgemask & MASK_HE_INNER;
    if ( !half0 || !half1 ) {
      switch ( cpdimen ) {
    case 2:
        for ( i = 0; i < nhe; i++ ) {
          CHECKIT1
          if ( ok ) {
            v0 = mhe[i].v0;  v1 = mhe[i].v1;
            MidPoint2d ( (point2d*)&mvc[2*v0], (point2d*)&mvc[2*v1],
                         (point2d*)midpoint );
            glVertex2dv ( midpoint );
            if ( half0 )
              glVertex2dv ( &mvc[2*v0] );
            else
              glVertex2dv ( &mvc[2*v1] );
          }
        }
        break;
    case 3:
        if ( spdimen == 2 ) {
          for ( i = 0; i < nhe; i++ ) {
            CHECKIT1
            if ( ok ) {
              v0 = mhe[i].v0;  v1 = mhe[i].v1;
              MidPoint3d ( (point3d*)&mvc[3*v0], (point3d*)&mvc[3*v1],
                           (point3d*)midpoint );
              glVertex4d ( midpoint[0], midpoint[1], 0.0, midpoint[2] );
              if ( half0 )
                glVertex4d ( mvc[3*v0], mvc[3*v0+1], 0.0, mvc[3*v0+2] );
              else
                glVertex4d ( mvc[3*v1], mvc[3*v1+1], 0.0, mvc[3*v1+2] );
            }
          }
        }
        else {
          for ( i = 0; i < nhe; i++ ) {
            CHECKIT1
            if ( ok ) {
              v0 = mhe[i].v0;  v1 = mhe[i].v1;
              MidPoint3d ( (point3d*)&mvc[3*v0], (point3d*)&mvc[3*v1],
                           (point3d*)midpoint );
              glVertex3dv ( midpoint );
              if ( half0 )
                glVertex3dv ( &mvc[3*v0] );
              else
                glVertex3dv ( &mvc[3*v1] );
            }
          }
        }
        break;
    case 4:
        for ( i = 0; i < nhe; i++ ) {
          CHECKIT1
          if ( ok ) {
            v0 = mhe[i].v0;  v1 = mhe[i].v1;
            MidPoint4d ( (point4d*)&mvc[4*v0], (point4d*)&mvc[4*v1],
                         (point4d*)midpoint );
            glVertex4dv ( midpoint );
            if ( half0 )
              glVertex4dv ( &mvc[4*v0] );
            else
              glVertex4dv ( &mvc[4*v1] );
          }
        }
        break;
    default:
        break;
      }
    }
    else {
      switch ( cpdimen ) {
    case 2:
        for ( i = 0; i < nhe; i++ ) {
          CHECKIT2
          if ( ok ) {
            glVertex2dv ( &mvc[2*mhe[i].v0] );
            glVertex2dv ( &mvc[2*mhe[i].v1] );
          }
        }
        break;
    case 3:
        if ( spdimen == 2 ) {
          for ( i = 0; i < nhe; i++ ) {
            CHECKIT2
            if ( ok ) {
              j = 3*mhe[i].v0;
              glVertex4d ( mvc[j], mvc[j+1], 0.0, mvc[j+2] );
              j = 3*mhe[i].v1;
              glVertex4d ( mvc[j], mvc[j+1], 0.0, mvc[j+2] );
            }
          }
        }
        else {
          for ( i = 0; i < nhe; i++ ) {
            CHECKIT2
            if ( ok ) {
              glVertex3dv ( &mvc[3*mhe[i].v0] );
              glVertex3dv ( &mvc[3*mhe[i].v1] );
            }
          }
        }
        break;
    case 4:
        for ( i = 0; i < nhe; i++ ) {
          CHECKIT2
          if ( ok ) {
            glVertex4dv ( &mvc[4*mhe[i].v0] );
            glVertex4dv ( &mvc[4*mhe[i].v1] );
          }
        }
        break;
    default:
        break;
      }
    }
  glEnd ();
#undef CHECKIT1
#undef CHECKIT2
} /*GeomObjectDrawMeshHalfedges*/

void GeomObjectDrawMeshVertex ( int cpdimen, int spdimen,
                                int nv, BSMvertex *mv, double *mvc, int *mvhei,
                                int nhe, BSMhalfedge *mhe,
                                int nfac, BSMfacet *mfac, int *mfhei,
                                int vertnum )
{
  glPointSize ( 5.0 );
  glBegin ( GL_POINTS );
  switch ( cpdimen ) {
case 2:
    glVertex2dv ( &mvc[2*vertnum] );
    break;
case 3:
    if ( spdimen == 3 )
      glVertex3dv ( &mvc[3*vertnum] );
    else
      glVertex4d ( mvc[3*vertnum], mvc[3*vertnum+1], 0.0, mvc[3*vertnum+2] );
    break;
case 4:
    glVertex4dv ( &mvc[4*vertnum] );
    break;
default:
    break;
  }
  glEnd ();
} /*GeomObjectDrawMeshVertex*/

