
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */ 
/* //////////////////////////////////////////////////// */

#define DATAGEN_DOM_PARAMS  3
#define DATAGEN_SURF_PARAMS 5

void InitGHKnotsf ( int hole_k, float *knots );
void InitGHVectors2f ( int hole_k, vector2f *v );
void InitGHVectors3f ( int hole_k, vector3f *v );
void LineInters2f ( const point2f *p0, const point2f *p1,
                    const point2f *q0, const point2f *q1,
                    point2f *p );
void LineInters3f ( const point3f *p0, const point3f *p1, 
                    const point3f *q0, const point3f *q1,
                    point3f *p );
int InitGHDomainNetf ( int hole_k, const vector2f *v,
                       int nparams, const float *param,
                       point2f *domcp );
int InitGHSurfNetf ( int hole_k, const vector3f *v,
                     int nparams, const float *param,
                     point3f *surfcp );

