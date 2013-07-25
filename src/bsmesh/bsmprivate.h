
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

int _bsm_AveragingFVN ( int inv, BSMvertex *imv, int *imvhei,
                        int inhe, BSMhalfedge *imhe,
                        int infac, BSMfacet *imfac, int *imfhei,
                        int fn, char *vtag );

boolean _bsm_RotateHalfedgei ( int degree, int *hei, int hn );

void _bsm_OutputWorkMeshd ( int spdimen,
                            int wnv, BSMvertex *wmv, int *wmvhei, double *wptc,
                            int wnhe, BSMhalfedge *wmhe,
                            int wnfac, BSMfacet *wmfac, int *wmfhei,
                            int *newvi, int *newhei, int *newfi,
                            int *onv, BSMvertex *omv, int *omvhei, double *optc,
                            int *onhe, BSMhalfedge *omhe,
                            int *onfac, BSMfacet *omfac, int *omfhei );

