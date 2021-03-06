
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>   
#include <stdio.h>
#include <math.h>  
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/gl.h>  
#include <GL/glu.h> 
#include <GL/glx.h>

#include "pkvaria.h" 
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "bsmesh.h"
#include "g2blendingd.h"
#include "egholed.h"
#include "mengerc.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

#include "editor.h"
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"

/* ////////////////////////////////////////////////////////////////////////// */
void GeomObjectBeginReading ( void *usrdata, int obj_type )
{
  rw_object_attributes *attrib;

  attrib = (rw_object_attributes*)usrdata;
        /* space dimension to be fixed later */
  switch ( obj_type ) {
case BSF_BEZIER_CURVE:
    attrib->go_being_read =
         (geom_object*)GeomObjectAddBezierCurve ( "", 2, false );
    goto init_attr;
case BSF_BEZIER_PATCH:
    attrib->go_being_read =
         (geom_object*)GeomObjectAddBezierPatch ( "", 2, false );
    goto init_attr;
case BSF_BSPLINE_CURVE:
    attrib->go_being_read =
         (geom_object*)GeomObjectAddBSplineCurve ( "", 2, false );
    goto init_attr;
case BSF_BSPLINE_PATCH:
    attrib->go_being_read =
         (geom_object*)GeomObjectAddBSplinePatch ( "", 2, false );
    goto init_attr;
case BSF_BSPLINE_MESH:
    attrib->go_being_read =
         (geom_object*)GeomObjectAddBSplineMesh ( "", 2, false );
    goto init_attr;
case BSF_BSPLINE_HOLE:  /* none processed yet */
    return;

      /* initialise the attributes to default values */
init_attr:
    attrib->obj_type = obj_type;
    attrib->integrity_ok = true;  /* in case of error this will be changed */
        /* points marking */
    attrib->nmk = 0;
    attrib->mk = NULL;
        /* halfedge marking */
    attrib->nhemk = 0;
    attrib->hemk = NULL;
        /* colour: */
    attrib->colour[0] = attrib->colour[1] = 1.0;
    attrib->colour[2] = 0.0;
        /* dependency identifiers */
    attrib->filedepname = -1;
    attrib->filedepnum = 0;
    attrib->filedepid = NULL;
        /* others - in future */
    break;

case BSF_TRIMMED_DOMAIN:
    attrib->ntrd = 0;
    attrib->trd_error = false;
    break;

default:
    break;
  }
} /*GeomObjectBeginReading*/

static void _GeomObjectGetCPMK ( geom_object *go, int *ncp, byte **mk )
{
  switch ( go->obj_type ) {
case GO_BEZIER_CURVE:
    *ncp = ((GO_BezierCurve*)go)->degree+1;
    *mk  = ((GO_BezierCurve*)go)->mkcp;
    break;
case GO_BEZIER_PATCH:
    *ncp = (((GO_BezierPatch*)go)->degree_u+1)*(((GO_BezierPatch*)go)->degree_v+1);
    *mk  = ((GO_BezierPatch*)go)->mkcp;
    break;
case GO_BSPLINE_CURVE:
    *ncp = ((GO_BSplineCurve*)go)->lastknot-((GO_BSplineCurve*)go)->degree;
    *mk  = ((GO_BSplineCurve*)go)->mkcp;
    break;
case GO_BSPLINE_PATCH:
    *ncp = (((GO_BSplinePatch*)go)->lastknot_u-((GO_BSplinePatch*)go)->degree_u)*
           (((GO_BSplinePatch*)go)->lastknot_v-((GO_BSplinePatch*)go)->degree_v);
    *mk  = ((GO_BSplinePatch*)go)->mkcp;
    break;
case GO_BSPLINE_MESH:
    *ncp = ((GO_BSplineMesh*)go)->nv;
    *mk  = ((GO_BSplineMesh*)go)->mkcp;
    break;
case GO_BSPLINE_HOLE:
    *ncp = 12*((GO_BSplineHole*)go)->hole_k + 1;
    *mk  = ((GO_BSplineHole*)go)->mkcp;
    break;
default:
    *ncp = 0;
    *mk = NULL;
    break;
  }
} /*_GeomObjectGetCPMK*/

static void _GeomObjectGetHEMK ( geom_object *go, int *nhe, byte **mk )
{
  if ( go->obj_type == GO_BSPLINE_MESH ) {
    *nhe = ((GO_BSplineMesh*)go)->nhe;
    *mk  = ((GO_BSplineMesh*)go)->mkhe;
  }
  else {
    *nhe = 0;
    *mk  = NULL;
  }
} /*_GeomObjectGetHEMK*/

void GeomObjectEndReading ( void *usrdata, int obj_type, boolean success )
{
  rw_object_attributes *attrib;
  byte                 *mk;
  int                  i, ncp, nhe, nelem;
  mbs_polycurved       *elem;

  attrib = (rw_object_attributes*)usrdata;
  if ( !attrib->go_being_read )
    return;
  switch ( obj_type ) {
case BSF_BEZIER_CURVE:
case BSF_BEZIER_PATCH:
case BSF_BSPLINE_CURVE:
case BSF_BSPLINE_PATCH:
case BSF_BSPLINE_MESH:
case BSF_BSPLINE_HOLE:
    if ( success ) {
      /* the object has been read in, it is the last in the list, */
      /* but if there is an error (e.g. mesh inconsistency), */
      /* the object has to be deleted from the list */
      if ( !attrib->integrity_ok ) {
        current_go = attrib->go_being_read;
        GeomObjectDeleteCurrent ();
        goto cleanup;
      }
      /* assign the attributes */
        /* point marking */
      _GeomObjectGetCPMK ( attrib->go_being_read, &ncp, &mk );
      if ( attrib->nmk && mk &&
           ncp > 0 && ncp == attrib->nmk )
        memcpy ( mk, attrib->mk, ncp*sizeof(byte) );
      else if ( ncp > 0 )
        memset ( mk, 0, ncp*sizeof(byte) );
      for ( i = 0; i < ncp; i++ )
        mk[i] |= MASK_CP_MOVEABLE;
      if ( obj_type == BSF_BSPLINE_MESH ) {
        /* halfedge marking */
        _GeomObjectGetHEMK ( attrib->go_being_read, &nhe, &mk );
        if ( attrib->nhemk && mk &&
             nhe > 0 && nhe == attrib->nhemk )
          memcpy ( mk, attrib->hemk, nhe*sizeof(byte) );
        else if ( nhe > 0 )
          memset ( mk, 0, nhe*sizeof(byte) );
      }
        /* colour */
      memcpy ( attrib->go_being_read->colour, attrib->colour, 3*sizeof(double) );
        /* move the dependencies */
      attrib->go_being_read->filedepname = attrib->filedepname;
      attrib->go_being_read->filedepnum = attrib->filedepnum;
      attrib->go_being_read->filedepid = attrib->filedepid;
        /* others - in future */
    }
    else {
      if ( attrib->filedepid ) free ( attrib->filedepid );
    }
cleanup:
    attrib->filedepname = -1;
    attrib->filedepnum = 0;
    attrib->filedepid = NULL;
    if ( attrib->mk ) {
      free ( attrib->mk );
      attrib->mk = NULL;
    }
    if ( attrib->hemk ) {
      free ( attrib->hemk );
      attrib->hemk = NULL;
    }
    attrib->go_being_read = NULL;
    break;

case BSF_TRIMMED_DOMAIN:
    if ( attrib->trdelem ) {
      if ( attrib->go_being_read &&
           !attrib->trd_error && attrib->ntrd > 0 ) {
        nelem = attrib->ntrd;
        elem = malloc ( nelem*sizeof(mbs_polycurved) );
        if ( elem ) {  /* pass the trimmed domain boundary to the object */
          memcpy ( elem, attrib->trdelem, nelem*sizeof(mbs_polycurved) );
          GeomObjectBSplinePatchEnterTrimmedDomain (
                   (GO_BSplinePatch*)attrib->go_being_read, nelem, elem );
          elem = attrib->trdelem;
        }
        else {  /* deallocate the trimmed domain boundary description */
          elem = attrib->trdelem;
          for ( i = 0; i < nelem; i++ ) {
            if ( elem[i].knots )  free ( elem[i].knots );
            if ( elem[i].points ) free ( elem[i].points );
          }
        }
        free ( elem );
      }
    }
    attrib->trdelem = NULL;
    attrib->ntrd = attrib->trdl = 0;
    break;

default:
    break;
  }
} /*GeomObjectEndReading*/

void GeomObjectReadColour ( void *usrdata, point3d *colour )
{
  rw_object_attributes *attrib;

  attrib = (rw_object_attributes*)usrdata;
  memcpy ( attrib->colour, colour, 3*sizeof(double) );
} /*GeomObjectReadColour*/

void GeomObjectReadCPMK ( void *usrdata, int ncp, unsigned int *mk )
{
  rw_object_attributes *attrib;
  int                  i;

  attrib = (rw_object_attributes*)usrdata;
  if ( !attrib->mk && ncp > 0 ) {
    attrib->mk = malloc ( ncp*sizeof(byte) );
    if ( attrib->mk ) {
      attrib->nmk = ncp;
      for ( i = 0; i < ncp; i++ )
        attrib->mk[i] = (byte)mk[i];
    }
  }
} /*GeomObjectReadCPMK*/

void GeomObjectReadHEMK ( void *usrdata, int nhe, unsigned int *mk )
{
  rw_object_attributes *attrib;
  int                  i;

  attrib = (rw_object_attributes*)usrdata;
  if ( !attrib->hemk && nhe > 0 ) {
    attrib->hemk = malloc ( nhe*sizeof(byte) );
    if ( attrib->hemk ) {
      attrib->nhemk = nhe;
      for ( i = 0; i < nhe; i++ )
        attrib->hemk[i] = (byte)mk[i];
    }
  }
} /*GeomObjectReadHEMK*/

boolean GeomObjectResolveDependencies ( void )
{
  geom_object *go;
  boolean     success;

  success = true;
  for ( go = first_go; go; go = go->next )
    if ( go->filedepid ) {
      switch ( go->obj_type ) {
    case GO_BEZIER_CURVE:
        success &= GeomObjectBCResolveDependencies ( (GO_BezierCurve*)go );
        break;
    case GO_BEZIER_PATCH:
        success &= GeomObjectBPResolveDependencies ( (GO_BezierPatch*)go );
        break;
    case GO_BSPLINE_CURVE:
        success &= GeomObjectBSCResolveDependencies ( (GO_BSplineCurve*)go );
        break;
    case GO_BSPLINE_PATCH:
        success &= GeomObjectBSPResolveDependencies ( (GO_BSplinePatch*)go );
        break;
    case GO_BSPLINE_MESH:
        success &= GeomObjectBSMResolveDependencies ( (GO_BSplineMesh*)go );
        break;
    case GO_BSPLINE_HOLE:
        success &= GeomObjectBSHResolveDependencies ( (GO_BSplineHole*)go );
        break;
    default:
        break;
      }
      free ( go->filedepid );
      go->filedepid = NULL;
      go->filedepnum = 0;
      go->filedepname = 0;
    }
  return success;
} /*GeomObjectResolveDependencies*/

void GeomObjectReadDependency ( void *usrdata,
                                int depname, int ndep, int *dep )
{
  rw_object_attributes *attrib;

  attrib = (rw_object_attributes*)usrdata;
  attrib->filedepname = depname;
  attrib->filedepnum = ndep;
  if ( ndep > 0 ) {
    attrib->filedepid = malloc ( ndep*sizeof(int) );
    if ( attrib->filedepid )
      memcpy ( attrib->filedepid, dep, ndep*sizeof(int) );
  }
  else
    attrib->filedepid = NULL;
} /*GeomObjectReadDependency*/

void GeomObjectReadTrimmedDomain ( void *usrdata, mbs_polycurved *elem )
{
  rw_object_attributes *attrib;
  mbs_polycurved       *trdelem;
  int                  trdl, ntrd;
  int                  dim, deg, lkn, ncp;

  attrib = (rw_object_attributes*)usrdata;
  if ( attrib->trd_error )
    return;
        /* reallocate, if necessary */
  ntrd = attrib->ntrd;
  if ( ntrd+1 >= attrib->trdl ) {
    trdl = attrib->ntrd+16;
    trdelem = malloc ( trdl*sizeof(mbs_polycurved) );
    if ( trdelem ) {
      if ( ntrd && attrib->trdelem ) {
        memcpy ( trdelem, attrib->trdelem, ntrd*sizeof(mbs_polycurved) );
        free ( attrib->trdelem );
      }
      attrib->trdelem = trdelem;
      attrib->trdl = trdl;
    }
    else {
      attrib->trd_error = true;
      return;
    }
  }
        /* copy the data from a buffer into the destination arrays */
  trdelem = attrib->trdelem;
  trdelem[ntrd].ident = elem->ident;
  trdelem[ntrd].closing = elem->closing;
  trdelem[ntrd].cpdimen = dim = elem->cpdimen;
  trdelem[ntrd].degree = deg = elem->degree;
  trdelem[ntrd].lastknot = lkn = elem->lastknot;
  trdelem[ntrd].points = NULL;
  if ( elem->knots ) {
    trdelem[ntrd].knots = malloc ( (lkn+1)*sizeof(double) );
    if ( trdelem[ntrd].knots ) {
      memcpy ( trdelem[ntrd].knots, elem->knots, (lkn+1)*sizeof(double) );
      ncp = lkn-deg;
    }
    else {
      attrib->trd_error = true;
      ncp = 0;
    }
  }
  else if ( deg > 1 ) {
    trdelem[ntrd].knots = NULL;
    ncp = deg+1;
  }
  else {
    trdelem[ntrd].knots = NULL;
    ncp = lkn+1;
  }
  if ( !attrib->trd_error ) {
    trdelem[ntrd].points = malloc ( ncp*dim*sizeof(double) );
    if ( trdelem[ntrd].points ) {
      memcpy ( trdelem[ntrd].points, elem->points, ncp*dim*sizeof(double) );
      attrib->ntrd = ntrd+1;
    }
    else
      attrib->trd_error = true;
  }
} /*GeomObjectReadTrimmedDomain*/

boolean GeomObjectReadFile ( char *filename, bsf_Camera_fptr CameraReader )
{
  bsf_UserReaders      readers;
  rw_object_attributes attrib;
  boolean              result;

        /* register procedures entering the objects read in */
  bsf_ClearReaders ( &readers );
  memset ( &attrib, 0, sizeof(rw_object_attributes) );
  readers.userData = (void*)&attrib;
  bsf_BeginReadingFuncd ( &readers, GeomObjectBeginReading );
  bsf_EndReadingFuncd ( &readers, GeomObjectEndReading );
  bsf_BC4ReadFuncd  ( &readers, GeomObjectReadBezierCurve, MAX_DEGREE );
  bsf_BSC4ReadFuncd ( &readers, GeomObjectReadBSplineCurve,
                      MAX_DEGREE, MAX_BSC_KNOTS-1 );
  bsf_BP4ReadFuncd  ( &readers, GeomObjectReadBezierPatch, MAX_DEGREE );
  bsf_BSP4ReadFuncd ( &readers, GeomObjectReadBSplinePatch,
                      MAX_DEGREE, MAX_BSP_KNOTS-1 );
  bsf_BSM4ReadFuncd ( &readers, GeomObjectReadBSplineMesh, MAX_DEGREE,
                      MAX_BSM_NV, MAX_BSM_NHE, MAX_BSM_NFAC );
  bsf_TrimmedReadFuncd ( &readers, GeomObjectReadTrimmedDomain );
/*
  bsf_BSH4ReadFuncd ( &readers, GeomObjectReadBSplineHole );
*/
  bsf_CPMarkReadFunc ( &readers, GeomObjectReadCPMK );
  bsf_HalfedgeMarkReadFunc ( &readers, GeomObjectReadHEMK );
  bsf_DependencyReadFunc ( &readers, GeomObjectReadDependency, MAX_DEPNUM );
  bsf_CameraReadFuncd ( &readers, CameraReader );
  bsf_ColourReadFuncd ( &readers, GeomObjectReadColour );
        /* read the file */
  result = bsf_ReadBSFiled ( filename, &readers );
  if ( attrib.mk )      /* in case of failure */
    free ( attrib.mk );
  result &= GeomObjectResolveDependencies ();
  return result;
} /*GeomObjectReadFile*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean GeomObjectWriteAttributes ( void *usrdata )
{
  void         *sp;
  geom_object  *go;
  unsigned int *mkcp;
  byte         *mk;
  int          ncp, nhe, i;

  sp = pkv_GetScratchMemTop ();
  go = (geom_object*)usrdata;
  if ( !go )
    goto failure;
        /* write the object's dependencies information and other */
        /* object-specific attributes */
  switch ( go->obj_type ) {
case GO_BEZIER_CURVE:
    GeomObjectWriteBCAttributes ( (GO_BezierCurve*)go );
    break;
case GO_BEZIER_PATCH:
    GeomObjectWriteBPAttributes ( (GO_BezierPatch*)go );
    break;
case GO_BSPLINE_CURVE:
    GeomObjectWriteBSCAttributes ( (GO_BSplineCurve*)go );
    break;
case GO_BSPLINE_PATCH:
    GeomObjectWriteBSPAttributes ( (GO_BSplinePatch*)go );
    break;
case GO_BSPLINE_MESH:
    GeomObjectWriteBSMAttributes ( (GO_BSplineMesh*)go );
    break;
case GO_BSPLINE_HOLE:
    GeomObjectWriteBSHAttributes ( (GO_BSplineHole*)go );
    break;
default:  /* this should never happen, but ignore just in case */
    break;
  }
        /* write the object's colour */
  if ( !bsf_WriteColour ( (void*)go->colour ) )
    goto failure;
        /* write the point markings */
  _GeomObjectGetCPMK ( go, &ncp, &mk );
  if ( mk ) {
    if ( ncp <= 0 )
      goto failure;
    mkcp = (unsigned int*)pkv_GetScratchMemi ( ncp );
    if ( !mkcp )
      goto failure;
    for ( i = 0; i < ncp; i++ )
      mkcp[i] = (unsigned int)mk[i];
    bsf_WritePointsMK ( ncp, mkcp );
    pkv_SetScratchMemTop ( sp );
  }
        /* write the halfedge markings - for meshes only */
  _GeomObjectGetHEMK ( go, &nhe, &mk );
  if ( mk ) {
    if ( nhe <= 0 )
      goto failure;
    mkcp = (unsigned int*)pkv_GetScratchMemi ( nhe );
    if ( !mkcp )
      goto failure;
    for ( i = 0; i < nhe; i++ )
      mkcp[i] = (unsigned int)mk[i];
    bsf_WriteHEdgesMK ( nhe, mkcp );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectWriteAttributes*/

boolean GeomObjectWriteObj ( geom_object *go )
{
  switch ( go->obj_type ) {
case GO_BEZIER_CURVE:
    if ( !GeomObjectWriteBezierCurve ( (GO_BezierCurve*)go ) )
      return false;
    break;
case GO_BEZIER_PATCH:
    if ( !GeomObjectWriteBezierPatch ( (GO_BezierPatch*)go ) )
      return false;
    break;
case GO_BSPLINE_CURVE:
    if ( !GeomObjectWriteBSplineCurve ( (GO_BSplineCurve*)go ) )
      return false;
    break;
case GO_BSPLINE_PATCH:
    if ( !GeomObjectWriteBSplinePatch ( (GO_BSplinePatch*)go ) )
      return false;
    break;
case GO_BSPLINE_MESH:
    if ( !GeomObjectWriteBSplineMesh ( (GO_BSplineMesh*)go ) )
      return false;
    break;
case GO_BSPLINE_HOLE:
    if ( !GeomObjectWriteBSplineHole ( (GO_BSplineHole*)go ) )
      return false;
    break;
default:
    break;
  }
  return true;
} /*GeomObjectWriteObj*/

typedef struct {
    int n, nmax;
    int *id;
  } ident_table;

static boolean _GeomObjectReadIdent ( void *usrdata, int objtype, int ident )
{
  ident_table *idt;

  idt = (ident_table*)usrdata;
  if ( idt->n >= idt->nmax ) {
    idt->id = realloc ( idt->id, (idt->nmax+256)*sizeof(int) );
    if ( idt->id )
      idt->nmax += 256;
    else
      return false;
  }
  idt->id[idt->n] = ident;
  idt->n ++;
  return true;
} /*_GeomObjectReadIdent*/

boolean GeomObjectMakeIdentifiers ( char *filename, boolean append )
{
  void        *sp;
  ident_table idt;
  geom_object *go;
  int         ident;
  int         i, j;

  sp = pkv_GetScratchMemTop ();
  idt.n = idt.nmax = 0;
  idt.id = NULL;
  if ( append ) {
        /* read the identifiers already present in the file */
        /* in order to avoid generating them */
    if ( !bsf_ReadIdentifiers ( filename, (void*)&idt, _GeomObjectReadIdent ) )
      goto failure;
    if ( idt.n > 1 )
      if ( pkv_SortFast ( sizeof(int), ID_SIGNED_INT, sizeof(int), 0,
                       idt.n, (void*)idt.id ) != SORT_OK )
        goto failure;
  }
        /* now assign identifiers to objects present in the dependency lists */
  for ( go = first_go; go; go = go->next )
    go->ident = -1;
  ident = 0;
  j = 0;
  for ( go = first_go; go; go = go->next )
    if ( go->dependencies ) {
      for ( i = 0; i < go->maxdn; i++ )
        if ( go->dependencies[i] ) {
          if ( go->dependencies[i]->ident == -1 ) {
              /* generate the next identifier */
            while ( j < idt.n && ident == idt.id[j] ) {
              ident ++;
              j ++;
            }
            go->dependencies[i]->ident = ident;
            ident ++;
          }
        }
    }

  if ( idt.id )
    free ( idt.id );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( idt.id )
    free ( idt.id );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectMakeIdentifiers*/

boolean GeomObjectWriteFile ( char *filename, char whattowrite,
                              boolean (*writeotherdata)( void *usrdata ),
                              void *usrdata,
                              boolean append )
{
  geom_object *go;

  if ( !GeomObjectMakeIdentifiers ( filename, append ) )
    return false;
  if ( !bsf_OpenOutputFile ( filename, append ) )
    return false;

  switch ( whattowrite ) {
case GO_WRITE_CURRENT:
    if ( !GeomObjectWriteObj ( current_go ) )
      goto failure;
    break;

case GO_WRITE_ACTIVE:
    for ( go = first_go; go; go = go->next )
      if ( go->active || go == current_go ) {
        if ( !GeomObjectWriteObj ( go ) )
          goto failure;
      }
    break;

case GO_WRITE_ALL:
    for ( go = first_go; go; go = go->next )
      if ( !GeomObjectWriteObj ( go ) )
        goto failure;
    break;

default:
    goto failure;
  }

  if ( writeotherdata ) {
    if ( !writeotherdata ( usrdata ) )
      goto failure;
  }

  bsf_CloseOutputFile ();
  return true;

failure:
  bsf_CloseOutputFile ();
  return false;
} /*GeomObjectWriteFile*/

