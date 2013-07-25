
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

#define MYINFINITY 1.0e+308

/*#define GRDIV(a,b) ((1.0-TAU)*(a)+TAU*(b)) */
#define GRDIV(a,b) (exp((1.0-TAU)*log(a)+TAU*log(b)))

#define SQUAREGDS (9*3)   /* square gradient data size */
#define SQUAREHDS (9*5*9) /* square Hessian data size */
#define DEG3 (6)
#define DEG2 (4)
#define DEG (2)
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

#define g11   Gstar[0]
#define g12   Gstar[1]
#define g22   Gstar[2]

#define tb11  Bstar[0]
#define tb12  Bstar[1]
#define tb22  Bstar[2]

#define dg11  Dgstar[0]
#define dg12  Dgstar[1]
#define dg22  Dgstar[2]

#define dtb11  Dbstar[0]
#define dtb12  Dbstar[1]
#define dtb22  Dbstar[2]

#define Ni10   Ni[0]
#define Ni01   Ni[1]
#define Ni20   Ni[2]
#define Ni11   Ni[3]
#define Ni02   Ni[4]

#define Mij0102 Mij[0]
#define Mij0111 Mij[1]
#define Mij0120 Mij[2]
#define Mij1001 Mij[3]
#define Mij1002 Mij[4]
#define Mij1011 Mij[5]
#define Mij1020 Mij[6]

#define Mij0110 -Mij1001
#define Mij0210 -Mij1002
#define Mij0201 -Mij0102
#define Mij1110 -Mij1011
#define Mij1101 -Mij0111
#define Mij2010 -Mij1020
#define Mij2001 -Mij0120

#define MijLen 7
int _g1bl_Asym3MatIndex ( int i, int j, int k, boolean *neg );
int _g1bl_SetupHessian1Profile ( int lastknotu, int lastknotv, int *prof );
int _g1bl_SetupClosedHessian1Profile ( int lastknotu, int lastknotv, int *prof );

boolean _g1bl_TabBasisFuncd ( int nkn, double **knots, double **coeff,
                              double **bf, double **dbf,
                              double **ddbf);

int _g1bl_Asym3MatIndex ( int i, int j, int k, boolean *neg );
double *_g1bl_NijIndd ( int nkn, double *Nijtab,
                        int i0, int i1, int j0, int j1,          
                        int l0, int l1 );

void g1bl_TabNid ( int nkn, double *bf, double *dbf, double *ddbf,
                   double *Nitab );
void g1bl_TabNijd ( int nkn, double *bf, double *dbf, double *ddbf,
                    double *Nijtab );
			     
double *_g1bl_MijIndd ( int nkn, double *Mijtab,
                        int i0, int i1, int j0, int j1, int l0, int l1 );
			
void g1bl_TabMijd ( int nkn, double *bf, double *dbf, double *ddbf,
                    double *Mijtab );		     
		     
void _g1bl_UCompPDerd ( int nkn, double *Nitab,
                        int pitch, point3d *cp,
                        int fcpn, int i, int j, vector3d *pder );
void _g1bl_UCompGStard ( const vector3d *pder, double *Gstar );
void _g1bl_UCompDGStard ( int nkn, double *Nitab,
                          int lastknotu, int lastknotv,
                          int ip0, int ip1, int jp0, int jp1,
                          int isq, int jsq, int i, int j,
                          const vector3d *pder, double *DGstar );
void _g1bl_UCompDDGStard ( int nkn, double *Nijtab,
                           int lastknotu, int lastknotv,
                           int ip0, int ip1, int jp0, int jp1,
                           int isq, int jsq, int i, int j,
                           double *DDGstar );
void _g1bl_UCompBStard ( const vector3d *pder, double *Bstar );
void _g1bl_UCompDBStard ( int nkn, double *Nitab,
                          int lastknotu, int lastknotv,
                          int ip0, int ip1, int jp0, int jp1,
                          int isq, int jsq, int i, int j,
                          const vector3d *pder, double *DBstar );

void _g1bl_UCompDDBStard ( int nkn, double *Mijtab,
                           int lastknotu, int lastknotv,
                           int ip0, int ip1, int jp0, int jp1,
                           int isq, int jsq, int i, int j,
                           const vector3d *pder, double *DDBstar );
			   
void g1bl_UFuncSQd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv,
                     int pitch, point3d *cp,
                     double tC,
                     int isq, int jsq,
                     double *ftab );
		     
void g1bl_biharmFuncSQd(int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv,
                     int pitch, point3d *cp,
                     double tC,
                     int isq, int jsq,
                     double *ftab);

void g1bl_QFuncSQd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv,
                     int pitch, point3d *cp,
                     double tC,
                     int isq, int jsq,
                     double *ftab );
		     
void g1bl_UFuncGradSQd ( int nkn, const double *qcoeff, double *Nitab,
                         int lastknotu, int lastknotv,
                         int pitch, point3d *cp,
                         double tC,
                         int isq, int jsq, int ip0, int ip1, int jp0, int jp1,
                         double *ftab, double *gtab );
void g1bl_UFuncGradHessianSQd ( int nkn, const double *qcoeff, double *Nitab,
                                double *Nijtab,double *Mijtab,
                                int lastknotu, int lastknotv,
                                int pitch, point3d *cp, double tC,
                                int isq, int jsq, int ip0, int ip1, int jp0, int jp1,
                                double *ftab, double *gtab, double *htab);

boolean _g1bl_ComputeDeltaQd ( int n, const int *prof, double **hrows,
                               const double *grad, const double *dcoeff,
                               double *dq );
boolean _g1bl_ShiftDecompHessiand ( int neqs, int hsize, int *prof, double *hessian,
                                    double *Lhessian, double **Lhrows, double nu );
double _g1bl_AuxNuFuncd ( int nknots,
                          const double *qcoeff, double *Nitab,
                          int neqs, int hsize, int *prof, double *hessian,
                          double *Lhessian, double **Lhrows, double nu,
                          int lastknotu, int lastknotv, int pitch,
                          point3d *cp, point3d *acp,
                          double *grad, double *dcoeff,
                          double tC, double *ftab );

double _g1bl_ClosedAuxNuFuncd ( int nknots,
                          const double *qcoeff, double *Nitab,
                          int neqs, int hsize, int *prof, double *hessian,
                          double *Lhessian, double **Lhrows, double nu,
                          int lastknotu, int lastknotv, int pitch,
                          point3d *cp, point3d *acp,
                          double *grad, double *dcoeff,
                          double tC, double *ftab );

boolean _g1bl_LazyHessiand ( int lastknotu, int lastknotv,
                     int ni, int nsq, point3d *acp, point3d *hcp,
                     char *dirtypt, char *dirtysq, boolean *all );

boolean _g1bl_ClosedLazyHessiand ( int lastknotu, int lastknotv,
                     int ni, int nsq, point3d *acp, point3d *hcp,
                     char *dirtypt, char *dirtysq, boolean *all );

