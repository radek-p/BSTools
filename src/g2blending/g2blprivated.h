
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#define _DEBUG

#define MYINFINITY 1.0e+308

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
#define puuu  pder[5]
#define puuv  pder[6] 
#define puvv  pder[7] 
#define pvvv  pder[8] 
#define g11   Gstar[0]
#define g12   Gstar[1]
#define g22   Gstar[2]
#define g11u  Gstar[3]  
#define g12u  Gstar[4]
#define g22u  Gstar[5]
#define g11v  Gstar[6]
#define g12v  Gstar[7]  
#define g22v  Gstar[8]
#define tb11  Bstar[0]
#define tb12  Bstar[1]
#define tb22  Bstar[2]
#define tb11u Bstar[3]
#define tb12u Bstar[4]
#define tb22u Bstar[5]
#define tb11v Bstar[6]
#define tb12v Bstar[7]
#define tb22v Bstar[8]
#define dg11  Dgstar[0]
#define dg12  Dgstar[1]
#define dg22  Dgstar[2]
#define dg11u Dgstar[3]
#define dg12u Dgstar[4]
#define dg22u Dgstar[5]
#define dg11v Dgstar[6]
#define dg12v Dgstar[7]
#define dg22v Dgstar[8]
#define dtb11  Dbstar[0]
#define dtb12  Dbstar[1]
#define dtb22  Dbstar[2]
#define dtb11u Dbstar[3]
#define dtb12u Dbstar[4]
#define dtb22u Dbstar[5]
#define dtb11v Dbstar[6]
#define dtb12v Dbstar[7]
#define dtb22v Dbstar[8]
#define Ni10   Ni[0]
#define Ni01   Ni[1]
#define Ni20   Ni[2]
#define Ni11   Ni[3]
#define Ni02   Ni[4]
#define Ni30   Ni[5]
#define Ni21   Ni[6]
#define Ni12   Ni[7]
#define Ni03   Ni[8]
#define Mij0102 Mij[0]
#define Mij0103 Mij[1]
#define Mij0111 Mij[2]
#define Mij0112 Mij[3]
#define Mij0120 Mij[4]
#define Mij0121 Mij[5]
#define Mij0130 Mij[6]
#define Mij1001 Mij[7]
#define Mij1002 Mij[8]
#define Mij1003 Mij[9]
#define Mij1011 Mij[10]
#define Mij1012 Mij[11]
#define Mij1020 Mij[12]
#define Mij1021 Mij[13]
#define Mij1030 Mij[14]
#define Mij1102 Mij[15]
#define Mij2002 Mij[16]
#define Mij2011 Mij[17]
#define Mij0201 -Mij0102
#define Mij0301 -Mij0103
#define Mij1101 -Mij0111
#define Mij1201 -Mij0112
#define Mij2001 -Mij0120
#define Mij2101 -Mij0121
#define Mij3001 -Mij0130
#define Mij0110 -Mij1001
#define Mij0210 -Mij1002
#define Mij0310 -Mij1003
#define Mij1110 -Mij1011
#define Mij1210 -Mij1012
#define Mij2010 -Mij1020
#define Mij2110 -Mij1021
#define Mij3010 -Mij1030
#define Mij0211 -Mij1102
#define Mij0220 -Mij2002
#define Mij1120 -Mij2011

int _g2bl_SetupHessian1Profile ( int lastknotu, int lastknotv, int *prof );
int _g2bl_SetupClosedHessian1Profile ( int lastknotu, int lastknotv, int *prof );

boolean _g2bl_TabBasisFuncd ( int nkn, double **knots, double **coeff,
                              double **bf, double **dbf,
                              double **ddbf, double **dddbf );

double *_g2bl_NijIndd ( int nkn, double *Nijtab,
                        int i0, int i1, int j0, int j1,          
                        int l0, int l1 );
double *_g2bl_MijIndd ( int nkn, double *Mijtab,
                        int i0, int i1, int j0, int j1, int l0, int l1 );

void g2bl_TabNid ( int nkn, double *bf, double *dbf, double *ddbf, double *dddbf,
                   double *Nitab );
void g2bl_TabNijd ( int nkn, double *bf, double *dbf, double *ddbf,
                    double *Nijtab );
void g2bl_TabMijd ( int nkn, double *bf, double *dbf, double *ddbf, double *dddbf,
                    double *Mijtab );

void _g2bl_UCompPDerd ( int nkn, double *Nitab,
                        int pitch, point3d *cp,
                        int fcpn, int i, int j, vector3d *pder );
void _g2bl_UCompGStard ( const vector3d *pder, double *Gstar );
void _g2bl_UCompDGStard ( int nkn, double *Nitab,
                          int lastknotu, int lastknotv,
                          int ip0, int ip1, int jp0, int jp1,
                          int isq, int jsq, int i, int j,
                          const vector3d *pder, double *DGstar );
void _g2bl_UCompDDGStard ( int nkn, double *Nijtab,
                           int lastknotu, int lastknotv,
                           int ip0, int ip1, int jp0, int jp1,
                           int isq, int jsq, int i, int j,
                           double *DDGstar );
void _g2bl_UCompBStard ( const vector3d *pder, double *Bstar );
void _g2bl_UCompDBStard ( int nkn, double *Nitab,
                          int lastknotu, int lastknotv,
                          int ip0, int ip1, int jp0, int jp1,
                          int isq, int jsq, int i, int j,
                          const vector3d *pder, double *DBstar );
void _g2bl_UCompDDBStard ( int nkn, double *Mijtab,
                           int lastknotu, int lastknotv,
                           int ip0, int ip1, int jp0, int jp1,
                           int isq, int jsq, int i, int j,
                           const vector3d *pder, double *DDBstar );

void _g2bl_UFuncSQIntegrandd ( vector3d pder[11], double *first, double *second );
void g2bl_UFuncSQd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv,
                     int pitch, point3d *cp,
                     double tC,
                     int isq, int jsq,
                     double *ftab );
void g2bl_UFuncGradSQd ( int nkn, const double *qcoeff, double *Nitab,
                         int lastknotu, int lastknotv,
                         int pitch, point3d *cp,
                         double tC,
                         int isq, int jsq, int ip0, int ip1, int jp0, int jp1,
                         double *ftab, double *gtab );
void g2bl_UFuncGradHessianSQd ( int nkn, const double *qcoeff, double *Nitab,
                                double *Nijtab, double *Mijtab,
                                int lastknotu, int lastknotv,
                                int pitch, point3d *cp,
                                double tC,
                                int isq, int jsq, int ip0, int ip1, int jp0, int jp1,
                                double *ftab, double *gtab, double *htab );

boolean _g2bl_ComputeDeltaQd ( int n, const int *prof, double **hrows,
                               const double *grad, const double *dcoeff,
                               double *dq );
boolean _g2bl_ShiftDecompHessiand ( int neqs, int hsize, int *prof, double *hessian,
                                    double *Lhessian, double **Lhrows, double nu );
double _g2bl_AuxNuFuncd ( int nknots,
                          const double *qcoeff, double *Nitab,
                          int neqs, int hsize, int *prof, double *hessian,
                          double *Lhessian, double **Lhrows, double nu,
                          int lastknotu, int lastknotv, int pitch,
                          point3d *cp, point3d *acp,
                          double *grad, double *dcoeff,
                          double tC, double *ftab );

double _g2bl_ClosedAuxNuFuncd ( int nknots,
                          const double *qcoeff, double *Nitab,
                          int neqs, int hsize, int *prof, double *hessian,
                          double *Lhessian, double **Lhrows, double nu,
                          int lastknotu, int lastknotv, int pitch,
                          point3d *cp, point3d *acp,
                          double *grad, double *dcoeff,
                          double tC, double *ftab );

boolean _g2bl_LazyHessiand ( int lastknotu, int lastknotv,
                     int ni, int nsq, point3d *acp, point3d *hcp,
                     char *dirtypt, char *dirtysq, boolean *all );

boolean _g2bl_ClosedLazyHessiand ( int lastknotu, int lastknotv,
                     int ni, int nsq, point3d *acp, point3d *hcp,
                     char *dirtypt, char *dirtysq, boolean *all );

boolean _g2mbl_DivideIntervald ( double *ga, double *gc, double *gd, double *gb,
                                 double *fa, double *fc, double *fd, double *fb );

