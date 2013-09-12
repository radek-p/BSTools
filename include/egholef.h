
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* this header file is intended to be #included by the #include files        */
/* eg1holef.h and eg2holef.h, which in turn may be #included by applications */

#ifndef EGHOLEF_H
#define EGHOLEF_H

#ifndef CONST_  /* a dirty trick to suppress many compiler warning messages */
#define CONST_ const
#endif

#ifdef __cplusplus   
extern "C" {
#endif

#define GH_MAX_K 16  /* maximal number of hole sides */

typedef struct GHoleDomainf {
    int     hole_k;         /* number of hole sides */
    float   *hole_knots;    /* knots of the surface with the hole */
    point2f *domain_cp;     /* control points of domain definition */
    boolean basisG1, basisG2;
    void    *privateG;      /* points to a common private data for G1 and G2 */
    void    *privateG1;     /* points to a G1HolePrivateRecf structure */
    void    *SprivateG1;    /* points to a G1HoleSPrivateRecf structire */
    void    *privateG2;     /* points to a G1HolePrivateRecf structure */
    void    *SprivateG2;    /* points to a G1HoleSPrivateRecf structire */
    void    *usrptr;        /* points to arbitrary application stuff */
    int     error_code;     /* explains possible failure reasons */
  } GHoleDomainf;


GHoleDomainf* gh_CreateDomainf ( int     hole_k,
                                 float   *hole_knots,
                                 point2f *domain_cp );
void gh_DestroyDomainf ( GHoleDomainf *domain );

boolean gh_GetBspInd ( int hole_k, int i, int j, int *ind );
int gh_DrawBFcpn ( int hole_k, unsigned char *bfcpn );

boolean gh_FindDomSurrndBezPatchesf ( GHoleDomainf *domain );
float gh_DomainDiamf ( GHoleDomainf *domain );
float gh_HoleDomainAreaf ( GHoleDomainf *domain, boolean symmetric );

/* ///////////////////////////////////////////////////////////////////////// */
/* drawing procedure */
void gh_DrawDomSurrndPatchesf ( GHoleDomainf *domain,
               void (*drawpatch) ( int n, int m, const point2f *cp ) );
void gh_GetDomSurrndBFuncf ( GHoleDomainf *domain, int fn, int i, int j,
                             float *bf );

/* ///////////////////////////////////////////////////////////////////////// */
/* domain net data, eigenvectors of the mesh refinement operators */
extern point2f egh_eigendom3f[], egh_eigendom5f[], egh_eigendom6f[],
               egh_eigendom7f[], egh_eigendom8f[], egh_eigendom9f[],
               egh_eigendom10f[], egh_eigendom11f[], egh_eigendom12f[],
               egh_eigendom13f[], egh_eigendom14f[], egh_eigendom15f[],
               egh_eigendom16f[], *egh_eigendomcpf[];
extern float egh_eigenvalf[GH_MAX_K-3];

#ifdef __cplusplus
}
#endif

#endif

