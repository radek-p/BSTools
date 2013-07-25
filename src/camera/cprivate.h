
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void _CameraUpdateRotCentref ( CameraRecf *CPos );
void _CameraUpdateRotCentred ( CameraRecd *CPos );
#ifdef STEREOF_H
void _StereoUpdateRotCentref ( StereoRecf *CPos );
#endif
#ifdef STEREOD_H
void _StereoUpdateRotCentred ( StereoRecd *CPos );
#endif

