
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

#ifndef STEREOD_H
#define STEREOD_H

#ifndef CAMERAD_H
#include "camerad.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct StereoRecd {
  point3d    position;         /* position of the observer        */
  double      d;               /* distance between the eyes       */
  double      l;               /* distance from the object centre */
  CameraRecd left, right;      /* cameras of the eyes             */
  trans3d    STr, STrInv;
} StereoRecd;


void StereoInitFramed ( StereoRecd *Stereo, boolean upside,
                        short width, short height, short xmin, short ymin,
                        double aspect, int ncplanes );
void StereoSetDimd ( StereoRecd *Stereo, double f, double d, double l );
void StereoSetMagd ( StereoRecd *Stereo, char mag );
void StereoSetDepthRanged ( StereoRecd *Stereo, double zmin, double zmax );
void StereoSetMappingd ( StereoRecd *Stereo );
void StereoInitPosd ( StereoRecd *Stereo );
void StereoSetRotCentred ( StereoRecd *Stereo,
                           point3d *centre,
                           boolean global_coord, boolean global_fixed );
void StereoUpdateRotCentred ( StereoRecd *Stereo );
void StereoMoveGd ( StereoRecd *Stereo, vector3d *v );
void StereoMoveCd ( StereoRecd *Stereo, vector3d *v );
void StereoRotGd ( StereoRecd *Stereo, double _psi, double _theta, double _phi );
#define StereoRotXGd(Stereo,angle) \
  StereoRotGd ( Stereo, 0.0, angle, 0.0 )
#define StereoRotYGd(Stereo,angle) \
  StereoRotGd ( Stereo, 0.5*PI, angle, -0.5*PI )
#define StereoRotZGd(Stereo,angle) \
  StereoRotGd ( Stereo, angle, 0.0, 0.0 )
void StereoRotVGd ( StereoRecd *Stereo, vector3d *v, double angle );
void StereoRotCd ( StereoRecd *Stereo, double _psi, double _theta, double _phi );
#define StereoRotXCd(Stereo,angle) \
  StereoRotCd ( Stereo, 0.0, angle, 0.0 )
#define StereoRotYCd(Stereo,angle) \
  StereoRotCd ( Stereo, 0.5*PI, angle, -0.5*PI )
#define StereoRotZCd(Stereo,angle) \
  StereoRotCd ( Stereo, angle, 0.0, 0.0 )
void StereoRotVCd ( StereoRecd *Stereo, vector3d *v, double angle );
void StereoZoomd ( StereoRecd *Stereo, double fchange );

#ifdef __cplusplus
}
#endif

#endif /*STEREOD_H*/

