
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* this file was written by Mateusz Markowski                                */
/* and modified by Przemyslaw Kiciak                                         */

#define MYINFINITY 1.0e+37

/*#define GRDIV(a,b) ((1.0-TAU)*(a)+TAU*(b)) */
#define GRDIV(a,b) (exp((1.0-TAU)*log(a)+TAU*log(b)))

#define SQUAREGDS (16*3)   /* square gradient data size */
#define SQUAREHDS (8*17*9) /* square Hessian data size */

/* squares may be considered "dirty", which means that the integrals */
/* related with the functional value, gradient or hessian need to be */
/* recalculated; below are bit masks for denoting that */
#define DIRTY_FUNC 0x01
#define DIRTY_GRAD 0x02
#define DIRTY_HESS 0x04

#define pu    pder[0]
#define pv    pder[1]
#define puu   pder[2]
#define puv   pder[3]
#define pvv   pder[4]
/*#define puuu  pder[5]
#define puuv  pder[6] 
#define puvv  pder[7] 
#define pvvv  pder[8] */
#define g11   Gstar[0]
#define g12   Gstar[1]
#define g22   Gstar[2]
/*#define g11u  Gstar[3]  
#define g12u  Gstar[4]
#define g22u  Gstar[5]
#define g11v  Gstar[6]
#define g12v  Gstar[7]  
#define g22v  Gstar[8]*/
#define tb11  Bstar[0]
#define tb12  Bstar[1]
#define tb22  Bstar[2]
/*#define tb11u Bstar[3]
#define tb12u Bstar[4]
#define tb22u Bstar[5]
#define tb11v Bstar[6]
#define tb12v Bstar[7]
#define tb22v Bstar[8]*/
#define dg11  Dgstar[0]
#define dg12  Dgstar[1]
#define dg22  Dgstar[2]
/*#define dg11u Dgstar[3]
#define dg12u Dgstar[4]
#define dg22u Dgstar[5]
#define dg11v Dgstar[6]
#define dg12v Dgstar[7]
#define dg22v Dgstar[8]*/
#define dtb11  Dbstar[0]
#define dtb12  Dbstar[1]
#define dtb22  Dbstar[2]
/*#define dtb11u Dbstar[3]
#define dtb12u Dbstar[4]
#define dtb22u Dbstar[5]
#define dtb11v Dbstar[6]
#define dtb12v Dbstar[7]
#define dtb22v Dbstar[8]*/
#define Ni10   Ni[0]
#define Ni01   Ni[1]
#define Ni20   Ni[2]
#define Ni11   Ni[3]
#define Ni02   Ni[4]
/*#define Ni30   Ni[5]
#define Ni21   Ni[6]
#define Ni12   Ni[7]
#define Ni03   Ni[8]*/

int _g2bl_Asym3MatIndex ( int i, int j, int k, boolean *neg );
int _g2bl_SetupHessian1Profile ( int lastknotu, int lastknotv, int *prof );

boolean _g2bl_TabBasisFuncf ( int nkn, float **knots, float **coeff,
                              float **bf, float **dbf,
                              float **ddbf, float **dddbf );

int _g2bl_Asym3MatIndex ( int i, int j, int k, boolean *neg );
float *_g2bl_NijIndf ( int nkn, float *Nijtab,
                        int i0, int i1, int j0, int j1,          
                        int l0, int l1 );
float *_g2bl_NijkIndf ( int nkn, float *Nijktab,
                         int i0, int i1, int j0, int j1, int k0, int k1,
                         int l0, int l1, boolean *neg );

void g2bl_TabNif ( int nkn, float *bf, float *dbf, float *ddbf, float *dddbf,
                   float *Nitab );
void g2bl_TabNijf ( int nkn, float *bf, float *dbf, float *ddbf,
                    float *Nijtab );
void g2bl_TabNijkf ( int nkn, float *bf, float *dbf, float *ddbf, float *dddbf,
                     float *Nijktab );

void _g2bl_UCompPDerf ( int nkn, float *Nitab,
                        int pitch, point3f *cp,
                        int fcpn, int i, int j, vector3f *pder );
void _g2bl_UCompGStarf ( const vector3f *pder, float *Gstar );
void _g2bl_UCompDGStarf ( int nkn, float *Nitab,
                          int lastknotu, int lastknotv,                
                          int isq, int jsq, int i, int j,       
                          const vector3f *pder, float *DGstar );
void _g2bl_UCompDDGStarf ( int nkn, float *Nijtab,
                           int lastknotu, int lastknotv,
                           int isq, int jsq, int i, int j,
                           float *DDGstar );
void _g2bl_UCompBStarf ( const vector3f *pder, float *Bstar );
void _g2bl_UCompDBStarf ( int nkn, float *Nitab,
                          int lastknotu, int lastknotv,
                          int isq, int jsq, int i, int j,
                          const vector3f *pder, float *DBstar );
void _g2bl_CompDDBStarf ( int nkn, float *Nijktab,
                          int lastknotu, int lastknotv,
                          int isq, int jsq, int i, int j,
                          int pitch, const point3f *cp,
                          float *DDBstar );

void g2bl_UFuncSQf ( int nkn, const float *qcoeff, float *Nitab,
                     int lastknotu, int lastknotv,
                     int pitch, point3f *cp,
                     float tC,
                     int isq, int jsq,
                     float *ftab );
void g2bl_UFuncGradSQf ( int nkn, const float *qcoeff, float *Nitab,
                         int lastknotu, int lastknotv,
                         int pitch, point3f *cp,
                         float tC,
                         int isq, int jsq,
                         float *ftab, float *gtab );
void g2bl_UFuncGradHessianSQf ( int nkn, const float *qcoeff, float *Nitab,
                                float *Nijtab, float *Nijktab,
                                int lastknotu, int lastknotv,
                                int pitch, point3f *cp,
                                float tC,
                                int isq, int jsq,
                                float *ftab, float *gtab, float *htab );

boolean _g2bl_ComputeDeltaQf ( int n, const int *prof, float **hrows,
                               const float *grad, const float *dcoeff,
                               float *dq );
boolean _g2bl_ShiftDecompHessianf ( int neqs, int hsize, int *prof, float *hessian,
                                    float *Lhessian, float **Lhrows, float nu );
float _g2bl_AuxNuFuncf ( int nknots,
                          const float *qcoeff, float *Nitab,
                          int neqs, int hsize, int *prof, float *hessian,
                          float *Lhessian, float **Lhrows, float nu,
                          int lastknotu, int lastknotv, int pitch,
                          point3f *cp, point3f *acp,
                          float *grad, float *dcoeff,
                          float tC, float *ftab );

boolean _g1bl_LazyHessianf ( int lastknotu, int lastknotv,
                     int ni, int nsq, point3f *acp, point3f *hcp,
                     char *dirtypt, char *dirtysq, boolean *all );

