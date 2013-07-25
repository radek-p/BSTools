
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2010                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* this header file is intended to be #included by the #include files        */
/* eg1holed.h and eg2holed.h, which in turn may be #included by applications */

#ifndef EGHOLED_H
#define EGHOLED_H

#ifndef CONST_  /* a dirty trick to suppress many compiler warning messages */
#define CONST_ const
#endif

#ifdef __cplusplus   
extern "C" {
#endif

#define GH_MAX_K 16  /* maximal number of hole sides */

typedef struct GHoleDomaind {
    int     hole_k;         /* number of hole sides */
    double   *hole_knots;    /* knots of the surface with the hole */
    point2d *domain_cp;     /* control points of domain definition */
    boolean basisG1, basisG2;
    void    *privateG;      /* points to a common private data dor G1 and G2 */
    void    *privateG1;     /* points to a G1HolePrivateRecd structure */
    void    *SprivateG1;    /* points to a G1HoleSPrivateRecd structire */
    void    *privateG2;     /* points to a G1HolePrivateRecd structure */
    void    *SprivateG2;    /* points to a G1HoleSPrivateRecd structire */
    void    *usrptr;        /* points to arbitrary application stuff */
    int     error_code;     /* explains possible failure reasons */
  } GHoleDomaind;


GHoleDomaind* gh_CreateDomaind ( int     hole_k,
                                 double   *hole_knots,
                                 point2d *domain_cp );
void gh_DestroyDomaind ( GHoleDomaind *domain );

void gh_GetBspInd ( int hole_k, int i, int j, int *ind );
int gh_DrawBFcpn ( int hole_k, unsigned char *bfcpn );

boolean gh_FindDomSurrndBezPatchesd ( GHoleDomaind *domain );
double gh_DomainDiamd ( GHoleDomaind *domain );
double gh_HoleDomainAread ( GHoleDomaind *domain, boolean symmetric );

/* ///////////////////////////////////////////////////////////////////////// */
/* drawing procedure */
void gh_DrawDomSurrndPatchesd ( GHoleDomaind *domain,
               void (*drawpatch) ( int n, int m, const point2d *cp ) );
void gh_GetDomSurrndBFuncd ( GHoleDomaind *domain, int fn, int i, int j,
                             double *bf );

/* ///////////////////////////////////////////////////////////////////////// */
/* domain net data, eigenvectors of the mesh refinement operators */
extern point2d egh_eigendom3d[], egh_eigendom5d[], egh_eigendom6d[],
               egh_eigendom7d[], egh_eigendom8d[], egh_eigendom9d[],
               egh_eigendom10d[], egh_eigendom11d[], egh_eigendom12d[],
               egh_eigendom13d[], egh_eigendom14d[], egh_eigendom15d[],
               egh_eigendom16d[], *egh_eigendomcpd[];
extern double egh_eigenvald[GH_MAX_K-3];

#ifdef __cplusplus
}
#endif

#endif

