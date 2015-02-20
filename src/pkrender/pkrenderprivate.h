
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

boolean BuildObjectTree ( pkRenderer *rend );
void DestroyTNodeTree ( rendertnode *node );

double ShapeFunc ( pkRenderer *rend, int func, ray3d *ray, renderobj *obj,
                   RayObjectIntersd *roint );
double cShapeFunc ( pkRenderer *rend, ray3d *ray, renderobj *obj,
                    RayObjectIntersd *roint );
boolean dShapeFunc ( pkRenderer *rend, ray3d *ray, renderobj *obj,
                     RayObjectIntersd *roint );

boolean FindRayTrInters ( pkRenderer *rend,
                          ray3d *ray, RayObjectIntersd *inters, int *k );
boolean IsInShadow ( pkRenderer *rend, point3d *p, vector3d *light );

void FindMinMaxShapeFunc ( pkRenderer *rend, int dens );

void GetTexColour ( pkRenderer *rend,
                    ray3d *ray, renderobj *obj, RayObjectIntersd *roint,
                    vector3d *texcolour );
void PhongLight ( pkRenderer *rend,
                  vector3d *tex, point3d *p, vector3d *nv,
                  vector3d *vv, vector3d *light, double light_intens,
                  vector3d *rgb );
void GetPixelColour ( pkRenderer *rend, ray3d *ray,
                      unsigned int *r, unsigned int *g, unsigned int *b );
