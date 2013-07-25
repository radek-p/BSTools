
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */ 
/* //////////////////////////////////////////////////// */

#define DATAGEN_DOM_PARAMS  3
#define DATAGEN_SURF_PARAMS 5

void InitGHKnotsd ( int hole_k, double *knots );
void InitGHVectors2d ( int hole_k, vector2d *v );
void InitGHVectors3d ( int hole_k, vector3d *v );
void LineInters2d ( const point2d *p0, const point2d *p1,
                    const point2d *q0, const point2d *q1,
                    point2d *p );
void LineInters3d ( const point3d *p0, const point3d *p1, 
                    const point3d *q0, const point3d *q1,
                    point3d *p );
int InitGHDomainNetd ( int hole_k, const vector2d *v,
                       int nparams, const double *param,
                       point2d *domcp );
int InitGHSurfNetd ( int hole_k, const vector3d *v,
                     int nparams, const double *param,
                     point3d *surfcp );

