
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#ifndef STEREOF_H
#define STEREOF_H

#ifndef CAMERAF_H
#include "cameraf.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct StereoRecf {
  point3f    position;         /* position of the observer        */
  float      d;                /* distance between the eyes       */
  float      l;                /* distance from the object centre */
  CameraRecf left, right;      /* cameras of the eyes             */
  trans3f    STr, STrInv;
} StereoRecf;


void StereoInitFramef ( StereoRecf *Stereo, boolean upside,
                        short width, short height, short xmin, short ymin,
                        float aspect, int ncplanes );
void StereoSetDimf ( StereoRecf *Stereo, float f, float d, float l );
void StereoSetMagf ( StereoRecf *Stereo, char mag );
void StereoSetDepthRangef ( StereoRecf *Stereo, float zmin, float zmax );
void StereoSetMappingf ( StereoRecf *Stereo );
void StereoInitPosf ( StereoRecf *Stereo );
void StereoSetRotCentref ( StereoRecf *Stereo,
                           point3f *centre,
                           boolean global_coord, boolean global_fixed );
void StereoUpdateRotCentref ( StereoRecf *Stereo );
void StereoMoveGf ( StereoRecf *Stereo, vector3f *v );
void StereoMoveCf ( StereoRecf *Stereo, vector3f *v );
void StereoRotGf ( StereoRecf *Stereo, float _psi, float _theta, float _phi );
#define StereoRotXGf(Stereo,angle) \
  StereoRotGf ( Stereo, 0.0, angle, 0.0 )
#define StereoRotYGf(Stereo,angle) \
  StereoRotGf ( Stereo, 0.5*PI, angle, -0.5*PI )
#define StereoRotZGf(Stereo,angle) \
  StereoRotGf ( Stereo, angle, 0.0, 0.0 )
void StereoRotVGf ( StereoRecf *Stereo, vector3f *v, float angle );
void StereoRotCf ( StereoRecf *Stereo, float _psi, float _theta, float _phi );
#define StereoRotXCf(Stereo,angle) \
  StereoRotCf ( Stereo, 0.0, angle, 0.0 )
#define StereoRotYCf(Stereo,angle) \
  StereoRotCf ( Stereo, 0.5*PI, angle, -0.5*PI )
#define StereoRotZCf(Stereo,angle) \
  StereoRotCf ( Stereo, angle, 0.0, 0.0 )
void StereoRotVCf ( StereoRecf *Stereo, vector3f *v, float angle );
void StereoZoomf ( StereoRecf *Stereo, float fchange );

#ifdef __cplusplus
}
#endif

#endif /*STEREOF_H*/

