
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MENUMASK        0xFFFF0000
#define MENUINC         0x00010000

/* window 0 widgets */
  /* top menu */
#define TOPMENU0                                 MENUINC
#define btnM00FILE                          (TOPMENU0+1)
#define btnM00EDIT                          (TOPMENU0+2)
#define btnM00TRANSFORM                     (TOPMENU0+3)
#define btnM00PICTURE                       (TOPMENU0+4)
#define btnM00ABOUT                         (TOPMENU0+5)

  /* side menu - edit */
#define SIDEMENU0                              2*MENUINC
#define swM01MKBIT_0                      (SIDEMENU0+ 1)
#define swM01MKBIT_1                      (SIDEMENU0+ 2)
#define swM01MKBIT_2                      (SIDEMENU0+ 3)
#define swM01MKBIT_3                      (SIDEMENU0+ 4)
#define swM01MKBIT_4                      (SIDEMENU0+ 5)
#define swM01MARK                         (SIDEMENU0+ 6)
#define swM01TRANSLATE                    (SIDEMENU0+ 7)
#define swM01SCALE                        (SIDEMENU0+ 8)
#define swM01ROTATE                       (SIDEMENU0+ 9)
#define swM01SHEAR                        (SIDEMENU0+10)
#define swM01PANZOOM                      (SIDEMENU0+11)
#define swM01SELECT_VERTEX                (SIDEMENU0+12)
#define swM01SELECT_EDGE                  (SIDEMENU0+13)
#define swM01COORDINATES                  (SIDEMENU0+14)
#define swM01STATUS                       (SIDEMENU0+15)
#define swM01COMMAND                      (SIDEMENU0+16)

  /* side menu - project */
#define SIDEMENU0e                             3*MENUINC
#define scwM01CONTENTS                   (SIDEMENU0e+ 1)
#define scwM01SCROLL                     (SIDEMENU0e+ 2)
#define textedM01ePX                     (SIDEMENU0e+ 3)
#define textedM01ePY                     (SIDEMENU0e+ 4)
#define textedM01ePZ                     (SIDEMENU0e+ 5)
#define textedM01eVX                     (SIDEMENU0e+ 6)
#define textedM01eVY                     (SIDEMENU0e+ 7)
#define textedM01eVZ                     (SIDEMENU0e+ 8)
#define textedM01eR                      (SIDEMENU0e+ 9)
#define textedM01eANGLE                  (SIDEMENU0e+10)
#define btnM01eTRANSLATE                 (SIDEMENU0e+11)
#define btnM01eSCALE                     (SIDEMENU0e+12)
#define btnM01eROTATE                    (SIDEMENU0e+13)
#define btnM01ePROJ_LINE                 (SIDEMENU0e+14)
#define btnM01ePROJ_PLANE                (SIDEMENU0e+15)
#define btnM01ePROJ_PLANE_UP             (SIDEMENU0e+16)
#define btnM01ePROJ_PLANE_DOWN           (SIDEMENU0e+17)
#define btnM01ePROJ_SPHERE               (SIDEMENU0e+18)
#define btnM01ePROJ_SPHERE_UP            (SIDEMENU0e+19)
#define btnM01ePROJ_SPHERE_DOWN          (SIDEMENU0e+20)
#define btnM01ePROJ_CYLINDER             (SIDEMENU0e+21)
#define btnM01ePROJ_CYLINDER_UP          (SIDEMENU0e+22)
#define btnM01ePROJ_CYLINDER_DOWN        (SIDEMENU0e+23)

  /* side menu - picture */
#define SIDEMENU0b                             4*MENUINC
#define btnM01bRENDER_INTERRUPT          (SIDEMENU0b+ 1)
#define swM01bPICTURE_ON                 (SIDEMENU0b+ 2)
#define btnM01bLIGHT                     (SIDEMENU0b+ 3)
#define swM01bGAUSSIAN_C                 (SIDEMENU0b+ 4)
#define swM01bMEAN_C                     (SIDEMENU0b+ 5)
#define swM01bLAMBERTISO_C               (SIDEMENU0b+ 6)
#define swM01bREFLECTION_C               (SIDEMENU0b+ 7)
#define swM01bHIGHLIGHT_C                (SIDEMENU0b+ 8)
#define swM01bSECTIONS_C                 (SIDEMENU0b+ 9)
#define swM01bGAUSSIAN_D                 (SIDEMENU0b+10)
#define swM01bMEAN_D                     (SIDEMENU0b+11)
#define swM01bLAMBERTISO_D               (SIDEMENU0b+12)
#define swM01bREFLECTION_D               (SIDEMENU0b+13)
#define swM01bHIGHLIGHT_D                (SIDEMENU0b+14)
#define swM01bSECTIONS_D                 (SIDEMENU0b+15)
#define slM01bREND_CFRANGE               (SIDEMENU0b+16)
#define slM01bREND_DFSF                  (SIDEMENU0b+17)
#define swM10bREFLECTION_FRAME           (SIDEMENU0b+18)
#define swM01bHIGHLIGHT_FRAME            (SIDEMENU0b+19)
#define swM01bSECTIONS_FRAME             (SIDEMENU0b+20)
#define swM01bANTIALIAS                  (SIDEMENU0b+21)
#define intwM01bRENDERINGTHREADS         (SIDEMENU0b+22)
  /* side menu - light */
#define btnM01bSHAPE_FUNCTION            (SIDEMENU0b+23)
#define swM01bLIGHT0DIR                  (SIDEMENU0b+24)
#define slM01bLIGHT0INT                  (SIDEMENU0b+25)
#define swM01bLIGHT1DIR                  (SIDEMENU0b+26)
#define slM01bLIGHT1INT                  (SIDEMENU0b+27)
#define swM01bLIGHT2DIR                  (SIDEMENU0b+28)
#define slM01bLIGHT2INT                  (SIDEMENU0b+29)
#define swM01bLIGHT3DIR                  (SIDEMENU0b+30)
#define slM01bLIGHT3INT                  (SIDEMENU0b+31)
#define slM01bLIGHTAMBIENTINT            (SIDEMENU0b+32)
#define swM01bSHADOWS                    (SIDEMENU0b+33)
#define btnM01bCAMERA                    (SIDEMENU0b+34)
  /* side menu - camera */
#define intwM01b_POS_X                   (SIDEMENU0b+35)
#define intwM01b_POS_Y                   (SIDEMENU0b+36)
#define intwM01b_POS_Z                   (SIDEMENU0b+37)
#define dialM01b_PSI                     (SIDEMENU0b+38)
#define dialM01b_THETA                   (SIDEMENU0b+39)
#define dialM01b_PHI                     (SIDEMENU0b+40)
#define intwM01b_PSI                     (SIDEMENU0b+41)
#define intwM01b_THETA                   (SIDEMENU0b+42)
#define intwM01b_PHI                     (SIDEMENU0b+43)
#define intwM01b_F                       (SIDEMENU0b+44)

  /* bottom menu */
#define BOTTOMMENU0                            5*MENUINC
#define textedM02COMMAND                 (BOTTOMMENU0+1)

  /* geometry window */
#define GEOMMENU0                              6*MENUINC

  /* File popup */
#define POPUP00                                7*MENUINC
#define btnP00OPEN                           (POPUP00+1)
#define btnP00SAVE                           (POPUP00+2)
#define btnP00SAVEAS                         (POPUP00+3)
#define btnP00EXIT                           (POPUP00+4)

  /* Open file popup */
#define POPUP01                                8*MENUINC
#define txtP01DIRSTR                         (POPUP01+1)
#define btnP01OPEN                           (POPUP01+2)
#define btnP01CANCEL                         (POPUP01+3)
#define lbP01DIRLIST                         (POPUP01+4)
#define lbP01FILELIST                        (POPUP01+5)
#define swP01HIDDEN                          (POPUP01+6)

  /* Save file as popup */
#define POPUP02                                9*MENUINC
#define txtP02OBJECTSTOSAVE                 0
#define txtP02DIRSTR                        (POPUP02+ 1)
#define btnP02SAVE                          (POPUP02+ 2)
#define btnP02CANCEL                        (POPUP02+ 3)
#define lbP02DIRLIST                        (POPUP02+ 4)
#define lbP02FILELIST                       (POPUP02+ 5)
#define txtedP02FILENAME                    (POPUP02+ 6)
#define swP02HIDDEN                         (POPUP02+ 7)
#define swP02ALL                            (POPUP02+ 8)
#define swP02ACTIVE                         (POPUP02+ 9)
#define swP02CURRENT                        (POPUP02+10)
#define swP02CAMERA                         (POPUP02+11)
#define swP02APPEND                         (POPUP02+12)

  /* Exit popup */
#define POPUP03                               10*MENUINC
#define btnP03EXIT                           (POPUP03+1)
#define btnP03SAVE                           (POPUP03+2)
#define btnP03CANCEL                         (POPUP03+3)

  /* 2D Geometry window */
#define GEOMWIN2D0                            12*MENUINC
  /* 3D Geometry window */
#define GEOMWIN3D0                            13*MENUINC


/* window 1 widgets */
  /* top menu */
#define TOPMENU1                              14*MENUINC
#define btnM10OBJECTS                       (TOPMENU1+1)
#define btnM10EDIT                          (TOPMENU1+2)
#define btnM10VIEW                          (TOPMENU1+3)
#define btnM10DATA                          (TOPMENU1+4)
#define btnM10OPTIONS                       (TOPMENU1+5)

  /* side menu */
#define SIDEMENU1                             15*MENUINC
#define swM11STATUS                        (SIDEMENU1+1)
#define swM11COMMAND                       (SIDEMENU1+2)

    /* widget sets specific for various geometric objects */
      /* Bezier curves */
#define SIDEMENU1_BEZC                        16*MENUINC
#define textedM1BEZC_NAME             (SIDEMENU1_BEZC+1)
#define intwM1BEZC_DEG                (SIDEMENU1_BEZC+2)

#define btnM1BEZC_VIEW_CURVE          (SIDEMENU1_BEZC+3)
#define btnM1BEZC_VIEW_CPOLY          (SIDEMENU1_BEZC+4)
#define swM1BEZC_VIEW_CURVATURE       (SIDEMENU1_BEZC+5)
#define slM1BEZC_SCALE_CURVATURE      (SIDEMENU1_BEZC+6)
#define swM1BEZC_VIEW_TORSION         (SIDEMENU1_BEZC+7)
#define slM1BEZC_SCALE_TORSION        (SIDEMENU1_BEZC+8)
#define intwM1BEZC_GRAPH_DENSITY      (SIDEMENU1_BEZC+9)
#define btnM1BEZC_COLOUR             (SIDEMENU1_BEZC+10)
#define slM1BEZC_PIPE_DIAMETER       (SIDEMENU1_BEZC+11)

      /* Bezier patches */
#define SIDEMENU1_BEZP                        17*MENUINC
#define textedM1BEZP_NAME             (SIDEMENU1_BEZP+1)
#define intwM1BEZP_DEGU               (SIDEMENU1_BEZP+2)
#define intwM1BEZP_DEGV               (SIDEMENU1_BEZP+3)
#define btnM1BEZP_FLIP                (SIDEMENU1_BEZP+4)

#define swM1BEZP_VIEW_SURF            (SIDEMENU1_BEZP+5)
#define swM1BEZP_VIEW_CNET            (SIDEMENU1_BEZP+6)
#define intwM1BEZP_DENSITY_U          (SIDEMENU1_BEZP+7)
#define intwM1BEZP_DENSITY_V          (SIDEMENU1_BEZP+8)
#define btnM1BEZP_COLOUR              (SIDEMENU1_BEZP+9)

      /* B-spline curves */
#define SIDEMENU1_BSC                         18*MENUINC
#define textedM1BSC_NAME               (SIDEMENU1_BSC+1)
#define intwM1BSC_DEG                  (SIDEMENU1_BSC+2)
#define btnM1BSC_UNIFORM               (SIDEMENU1_BSC+3)
#define btnM1BSC_REFINE                (SIDEMENU1_BSC+4)
#define swM1BSC_CLOSED                 (SIDEMENU1_BSC+5)
#define swM1BSC_MOVE_MANY_KNOTS        (SIDEMENU1_BSC+6)
#define swM1BSC_DOMAIN_COORD           (SIDEMENU1_BSC+7)
#define swM1BSC_DOMAIN_PANZOOM         (SIDEMENU1_BSC+8)

#define swM1BSC_VIEW_CURVE             (SIDEMENU1_BSC+9)
#define swM1BSC_VIEW_CPOLY            (SIDEMENU1_BSC+10)
#define swM1BSC_VIEW_BPOLY            (SIDEMENU1_BSC+11)
#define swM1BSC_VIEW_CURVATURE        (SIDEMENU1_BSC+12)
#define slM1BSC_SCALE_CURVATURE       (SIDEMENU1_BSC+13)
#define swM1BSC_VIEW_TORSION          (SIDEMENU1_BSC+14)
#define slM1BSC_SCALE_TORSION         (SIDEMENU1_BSC+15)
#define intwM1BSC_GRAPH_DENSITY       (SIDEMENU1_BSC+16)
#define btnM1BSC_COLOUR               (SIDEMENU1_BSC+17)
#define slM1BSC_PIPE_DIAMETER         (SIDEMENU1_BSC+18)

#define swM1MSC_MENGERC               (SIDEMENU1_BSC+19)
#define slM1MSC_MENGERC_EXP           (SIDEMENU1_BSC+20)
#define slM1MSC_MENGERC_P1            (SIDEMENU1_BSC+21) /* this and next 4 */
#define slM1MSC_MENGERC_P2            (SIDEMENU1_BSC+22) /* must be consecutive */
#define slM1MSC_MENGERC_P3            (SIDEMENU1_BSC+23)
#define slM1MSC_MENGERC_P4            (SIDEMENU1_BSC+24)
#define slM1MSC_MENGERC_P5            (SIDEMENU1_BSC+25)
#define intwM1BSC_MENGERC_QKN         (SIDEMENU1_BSC+26)
#define intwM1BSC_MENGERC_POPT        (SIDEMENU1_BSC+27)
#define intwM1BSC_MENGERC_MAXIT       (SIDEMENU1_BSC+28)
#define intwM1BSC_MENGERC_NPTHREADS   (SIDEMENU1_BSC+29)
#define swM1BSC_MENGERC_LOG           (SIDEMENU1_BSC+30)
#define btnM1BSC_MENGERC_OPTIMIZE     (SIDEMENU1_BSC+31)

      /* B-spline patches */
#define SIDEMENU1_BSP                         19*MENUINC
#define textedM1BSP_NAME              (SIDEMENU1_BSP+ 1)
#define btnM1BSP_TYPE                 (SIDEMENU1_BSP+ 2)
#define intwM1BSP_DEGU                (SIDEMENU1_BSP+ 3)
#define intwM1BSP_DEGV                (SIDEMENU1_BSP+ 4)
#define btnM1BSP_UNIFORM_U            (SIDEMENU1_BSP+ 5)
#define btnM1BSP_UNIFORM_V            (SIDEMENU1_BSP+ 6)
#define swM1BSP_CLOSED_U              (SIDEMENU1_BSP+ 7)
#define swM1BSP_CLOSED_V              (SIDEMENU1_BSP+ 8)
#define btnM1BSP_FLIP                 (SIDEMENU1_BSP+ 9)
#define swM1BSP_MOVE_MANY_KNOTS       (SIDEMENU1_BSP+10)
#define swM1BSP_DOMAIN_COORD          (SIDEMENU1_BSP+11)
#define swM1BSP_DOMAIN_PANZOOM        (SIDEMENU1_BSP+12)

#define swM1BSP_VIEW_SURF             (SIDEMENU1_BSP+13)
#define swM1BSP_VIEW_CNET             (SIDEMENU1_BSP+14)
#define intwM1BSP_DENSITY_U           (SIDEMENU1_BSP+15)
#define intwM1BSP_DENSITY_V           (SIDEMENU1_BSP+16)
#define btnM1BSP_COLOUR               (SIDEMENU1_BSP+17)

#define swM1BSP_BLENDING_CLAMPED      (SIDEMENU1_BSP+18)
#define swM1BSP_NHARMONIC             (SIDEMENU1_BSP+19)
#define btnM1BSP_PRETRANSFORMATION    (SIDEMENU1_BSP+20)
#define swM1BSP_BLENDING_ENTIRE       (SIDEMENU1_BSP+21)
#define intwM1BSP_G1BLENDING_UMIN     (SIDEMENU1_BSP+22)
#define intwM1BSP_G1BLENDING_UMAX     (SIDEMENU1_BSP+23)
#define intwM1BSP_G1BLENDING_VMIN     (SIDEMENU1_BSP+24)
#define intwM1BSP_G1BLENDING_VMAX     (SIDEMENU1_BSP+25)
#define intwM1BSP_G2BLENDING_UMIN     (SIDEMENU1_BSP+26)
#define intwM1BSP_G2BLENDING_UMAX     (SIDEMENU1_BSP+27)
#define intwM1BSP_G2BLENDING_VMIN     (SIDEMENU1_BSP+28)
#define intwM1BSP_G2BLENDING_VMAX     (SIDEMENU1_BSP+29)
#define slM1BSP_BLENDING_CPARAM       (SIDEMENU1_BSP+30)
#define intwM1BSP_NKN1                (SIDEMENU1_BSP+31)
#define intwM1BSP_NKN2                (SIDEMENU1_BSP+32)
#define intwM1BSP_MAXIT               (SIDEMENU1_BSP+33)
#define btnM1BSP_BLENDING_OPTIMIZE    (SIDEMENU1_BSP+34)
#define btnM1BSP_BLENDING_REFINE      (SIDEMENU1_BSP+35 )

#define swM1BSP_SPRODUCT_EQUATOR      (SIDEMENU1_BSP+36)
#define textedM1BSP_SPRODUCT_EQNAME   (SIDEMENU1_BSP+37)
#define intWM1BSP_SPRODUCT_EQDEG      (SIDEMENU1_BSP+38)
#define swM1BSP_SPRODUCT_EQRATIONAL   (SIDEMENU1_BSP+39)
#define swM1BSP_SPRODUCT_EQCLOSED     (SIDEMENU1_BSP+40)
#define swM1BSP_SPRODUCT_EQUNIFORM    (SIDEMENU1_BSP+41)
#define swM1BSP_SPRODUCT_MERIDIAN     (SIDEMENU1_BSP+42)
#define textedM1BSP_SPRODUCT_MERNAME  (SIDEMENU1_BSP+43)
#define intWM1BSP_SPRODUCT_MERDEG     (SIDEMENU1_BSP+44)
#define swM1BSP_SPRODUCT_MERRATIONAL  (SIDEMENU1_BSP+45)
#define swM1BSP_SPRODUCT_MERCLOSED    (SIDEMENU1_BSP+46)
#define swM1BSP_SPRODUCT_MERUNIFORM   (SIDEMENU1_BSP+47)

      /* B-spline meshes */
#define SIDEMENU1_BSM                         20*MENUINC
#define scwM1BSM_ECONTENTS            (SIDEMENU1_BSM+ 1)
#define scwM1BSM_ESCROLL              (SIDEMENU1_BSM+ 2)   
#define textedM1BSM_NAME              (SIDEMENU1_BSM+ 3)
#define intwM1BSM_DEG                 (SIDEMENU1_BSM+ 4)
#define btnM1BSM_REFINEMENT           (SIDEMENU1_BSM+ 5)
#define btnM1BSM_DOUBLING             (SIDEMENU1_BSM+ 6)
#define btnM1BSM_AVERAGING            (SIDEMENU1_BSM+ 7)
#define btnM1BSM_EXTRACTSUBMESH       (SIDEMENU1_BSM+ 8)

#define intwM1BSM_VERTEX0             (SIDEMENU1_BSM+ 9)
#define intwM1BSM_VERTEX1             (SIDEMENU1_BSM+10)
#define btnM1BSM_MARK_VERT            (SIDEMENU1_BSM+11)
#define btnM1BSM_UNMARK_VERT          (SIDEMENU1_BSM+12)
#define btnM1BSM_ENTER_LINE           (SIDEMENU1_BSM+13)
#define btnM1BSM_FILTER               (SIDEMENU1_BSM+14)
#define btnM1BSM_REMOVE_VERTEX        (SIDEMENU1_BSM+15)
#define btnM1BSM_DIVIDE_FACET         (SIDEMENU1_BSM+16)
#define intwM1BSM_EDGE0               (SIDEMENU1_BSM+17)
#define intwM1BSM_EDGE1               (SIDEMENU1_BSM+18)
#define btnM1BSM_SHRINK_EDGE          (SIDEMENU1_BSM+19)
#define btnM1BSM_CONTRACT_EDGE        (SIDEMENU1_BSM+20)
#define btnM1BSM_GLUE_EDGES           (SIDEMENU1_BSM+21)
#define btnM1BSM_GLUE_EDGE_LOOPS      (SIDEMENU1_BSM+22)
#define btnM1BSM_SEAL_HOLE            (SIDEMENU1_BSM+23)
#define btnM1BSM_SPLIT_BOUNDARY_EDGE  (SIDEMENU1_BSM+24)
#define intwM1BSM_FACET0              (SIDEMENU1_BSM+25)
#define intwM1BSM_FACET1              (SIDEMENU1_BSM+26)
#define btnM1BSM_REMOVE_FACET         (SIDEMENU1_BSM+27)
#define btnM1BSM_DOUBLE_FAC_EDGES     (SIDEMENU1_BSM+28)

#define swM1BSM_VIEW_SURF             (SIDEMENU1_BSM+29)
#define swM1BSM_VIEW_CNET             (SIDEMENU1_BSM+30)
#define intwM1BSM_DENSITY             (SIDEMENU1_BSM+31)
#define swM1BSM_VIEW_SPECIAL          (SIDEMENU1_BSM+32)
#define swM1BSM_HOLE_FILLING          (SIDEMENU1_BSM+33)
#define btnM1BSM_COLOUR               (SIDEMENU1_BSM+34)

#define intwM1BSM_DATA_KGON           (SIDEMENU1_BSM+35)
#define btnM1BSM_DATA_KGON            (SIDEMENU1_BSM+36)
#define btnM1BSM_DATA_TETRAHEDRON     (SIDEMENU1_BSM+37)
#define btnM1BSM_DATA_CUBE            (SIDEMENU1_BSM+38)
#define btnM1BSM_DATA_DODECAHEDRON    (SIDEMENU1_BSM+39)
#define intwM1BSM_DATA_KPRISM         (SIDEMENU1_BSM+40)
#define btnM1BSM_DATA_KPRISM          (SIDEMENU1_BSM+41)
#define swM1BSM_DATA_REPLACE          (SIDEMENU1_BSM+42)
#define swM1BSM_DATA_ADD              (SIDEMENU1_BSM+43)

#define scwM1BSM_OCONTENTS            (SIDEMENU1_BSM+44)
#define scwM1BSM_OSCROLL              (SIDEMENU1_BSM+45)   
#define swM1BSM_SUBDIVISION           (SIDEMENU1_BSM+46)
#define swM1BSM_BLENDING              (SIDEMENU1_BSM+47)
#define swM1BSM_G1                    (SIDEMENU1_BSM+48)
#define swM1BSM_G2                    (SIDEMENU1_BSM+49)
#define swM1BSM_G1_QUASI_G2           (SIDEMENU1_BSM+50)
#define slM1BSM_G1Q2_PARAM            (SIDEMENU1_BSM+51)
#define swM1BSM_HOLEFILL_COONS        (SIDEMENU1_BSM+52)
#define swM1BSM_HOLEFILL_BEZIER       (SIDEMENU1_BSM+53)
#define btnM1BSM_PRETRANSFORMATION    (SIDEMENU1_BSM+54)
#define slM1BSM_BLENDING_CPARAM       (SIDEMENU1_BSM+55)
#define swM1BSM_BLENDING_CONSTRAINTS  (SIDEMENU1_BSM+56)
#define intwM1BSM_NKN1                (SIDEMENU1_BSM+57)
#define intwM1BSM_NKN2                (SIDEMENU1_BSM+58)
#define intwM1BSM_MAXIT               (SIDEMENU1_BSM+59)
#define intwM1BSM_NLEVELS             (SIDEMENU1_BSM+60)
#define btnM1SUGGESTLEVELS            (SIDEMENU1_BSM+61)
#define btnM1SUGGESTBLOCKS            (SIDEMENU1_BSM+62)
#define intwM1BSM_NBLOCKS             (SIDEMENU1_BSM+63)
#define swM1BSM_SHAPE_ONLY            (SIDEMENU1_BSM+64)
#define swM1BSM_ALT_MULTILEVEL        (SIDEMENU1_BSM+65)
#define swM1BSM_COARSE_PRECOND        (SIDEMENU1_BSM+66)
#define textedM1BSM_COARSE_NAME       (SIDEMENU1_BSM+67)
#define intwM1BSM_NPTHREADS           (SIDEMENU1_BSM+68)
#define btnM1BSM_BLENDING_OPTIMIZE    (SIDEMENU1_BSM+69)
#define swM1BSM_LOG_IT                (SIDEMENU1_BSM+70)
#define btnM1BSM_OPTIMIZE_SPECIALS    (SIDEMENU1_BSM+71)

      /* B-spline holes */
#define SIDEMENU1_BSH                         21*MENUINC

#define btnM1BSH_VIEW_SURF             (SIDEMENU1_BSM+1)
#define btnM1BSH_VIEW_CNET             (SIDEMENU1_BSM+2)

  /* bottom menu */
#define BOTTOMMENU1                           22*MENUINC
#define textedM12COMMAND                 (BOTTOMMENU1+1)

  /* geometry menu */
#define GEOMMENU1                             23*MENUINC
#define GEOMWIN1_2D                        (GEOMMENU1+1)
#define GEOMWIN1_KN                        (GEOMMENU1+2)
#define GEOMWIN1_T2KN                      (GEOMMENU1+3)
#define GEOMWIN1_2DEQMER                   (GEOMMENU1+4)
#define GEOMWIN1_KNEQMER                   (GEOMMENU1+5)

  /* new object menu */
#define POPUP10                               24*MENUINC
#define btnP10BEZCURVE                       (POPUP10+1)
#define btnP10BEZPATCH                       (POPUP10+2)
#define btnP10BSCURVE                        (POPUP10+3)
#define btnP10BSPATCH                        (POPUP10+4)
#define btnP10BSMESH                         (POPUP10+5)
#define btnP10BSHOLE                         (POPUP10+6)

  /* New object popup; there are alternative widget sets */
  /* for various object types */
#define POPUP11                               25*MENUINC
#define txtedP11OBJNAME                     (POPUP11+ 1)
#define btnP11ADD_BEZC                      (POPUP11+ 2)
#define btnP11ADD_BEZP                      (POPUP11+ 3)
#define btnP11ADD_BSC                       (POPUP11+ 4)
#define btnP11ADD_BSP                       (POPUP11+ 5)
#define btnP11ADD_BSM                       (POPUP11+ 6)
#define btnP11ADD_BSH                       (POPUP11+ 7)
#define btnP11CANCEL                        (POPUP11+ 8)
#define swP11_2D                            (POPUP11+ 9)
#define swP11_3D                            (POPUP11+10)
#define swP11_RATIONAL                      (POPUP11+11)

  /* Object list popup */
#define POPUP12                               26*MENUINC
#define txtedP12OBJNAME                      (POPUP12+1)
#define lbP12OBJLIST                         (POPUP12+2)
#define btnP12NEW                            (POPUP12+3)
#define btnP12COPY                           (POPUP12+4)
#define btnP12DELETE                         (POPUP12+5)
#define btnP12PURGE                          (POPUP12+6)
#define swP12SHOW                            (POPUP12+7)
#define btnP12OK                             (POPUP12+8)

  /* Colour popup */
#define POPUP13                               27*MENUINC
#define spwP13_COLOUR                        (POPUP13+1)
#define btnP13_OK                            (POPUP13+2)
#define slP13_RED                            (POPUP13+4) /* % 4 == 0 */
#define slP13_GREEN                          (POPUP13+5) /* % 4 == 1 */
#define slP13_BLUE                           (POPUP13+6) /* % 4 == 2 */

  /* Pretransformation popup */
#define POPUP14                               28*MENUINC
#define txtedP14_A11                        (POPUP14+ 1) /* these 12 identifiers */
#define txtedP14_A12                        (POPUP14+ 2) /* must be consecutive */
#define txtedP14_A13                        (POPUP14+ 3)
#define txtedP14_A14                        (POPUP14+ 4)
#define txtedP14_A21                        (POPUP14+ 5)
#define txtedP14_A22                        (POPUP14+ 6)
#define txtedP14_A23                        (POPUP14+ 7)
#define txtedP14_A24                        (POPUP14+ 8)
#define txtedP14_A31                        (POPUP14+ 9)
#define txtedP14_A32                        (POPUP14+10)
#define txtedP14_A33                        (POPUP14+11)
#define txtedP14_A34                        (POPUP14+12)
#define btnP14_RESET                        (POPUP14+13)
#define btnP14_OK                           (POPUP14+14)

/* B-spline patch options popup */
#define POPUP15                               29*MENUINC
#define btnP15GENERAL                       (POPUP15+ 1)
#define btnP15SWEPT                         (POPUP15+ 2)
#define btnP15SPHERICAL                     (POPUP15+ 3)
#define btnP15LOFTED                        (POPUP15+ 4)
#define btnP15BLENDINGG1                    (POPUP15+ 5)
#define btnP15BLENDINGG2                    (POPUP15+ 6)


extern char txtNull[];
extern char txtFile[];
extern char txtObjects[];
extern char txtEdit[];
extern char txtTransform[];
extern char txtPicture[];
extern char txtAbout[];
extern char txtBit[];
extern char txtMarkUnmark[];
extern char txtMark[];
extern char txtUnmark[];
extern char txtTranslate[];
extern char txtScale[];
extern char txtRotate[];
extern char txtShear[];
extern char txtPanZoom[];
extern char txtCoordinates[];
extern char txtOK[];
extern char txtCancel[];
extern char txtExit[];
extern char txtHidden[];
extern char txtOpen[];
extern char txtSave[];
extern char txtSaveAs[];
extern char txtObjectsToSave[];
extern char txtAll[];
extern char txtCurrent[];
extern char txtAppend[];
extern char txtOptions[];
extern char txtView[];
extern char txtData[];
extern char txtName[];
extern char txtAdd[];
extern char txtReplace[];
extern char txtNew[];
extern char txtCopy[];
extern char txtDelete[];
extern char txtPurge[];
extern char txtActive[];
extern char txtBezierCurve[];
extern char txtBSplineCurve[];
extern char txtBezierPatch[];
extern char txtBSplinePatch[];
extern char txtBSplineMesh[];
extern char txtBSplineHole[];
extern char txtDegree[];
extern char txtDegreeU[];
extern char txtDegreeV[];
extern char txt2D[];
extern char txt3D[];
extern char txtRational[];
extern char txtRefine[];
extern char txtDouble[];
extern char txtAverage[];
extern char txtTetrahedron[];
extern char txtCube[];
extern char txtDodecahedron[];
extern char txt_gon[];
extern char txt_gonalPrism[];
extern char txtControlPolygon[];
extern char txtControlNet[];
extern char txtCurve[];
extern char txtBezierPolygons[];
extern char txtCurvatureGraph[];
extern char txtTorsionGraph[];
extern char txtGraphDensity[];
extern char txtSurface[];
extern char txtHoleFilling[];
extern char txtVertex[];
extern char txtLine[];
extern char txtEdge[];
extern char txtFacet[];
extern char txtRemove[];
extern char txtFilter[];
extern char txtDoubleEdges[];
extern char txtGlueEdges[];
extern char txtGlueLoops[];
extern char txtSealHole[];
extern char txtSplitEdge[];
extern char txtDivideFacet[];
extern char txtShrink[]; 
extern char txtContract[];
extern char txtDensity[];
extern char txtDensityU[];
extern char txtDensityV[];
extern char txtLevel[];
extern char txtSubdivision[];
extern char txtBlending[];
extern char txtSpecialNets[];
extern char txtUniform[];
extern char txtUniformU[];
extern char txtUniformV[];
extern char txtClosed[];
extern char txtClosedU[];
extern char txtClosedV[];
extern char txtFlip[];
extern char txtCoons[];
extern char txtBezier[];
extern char txtQKnots[];
extern char txtMaxIter[];
extern char txtStartFrom[];
extern char txtLevels[];
extern char txtBlocks[];
extern char txtConstraints[];
extern char txtShapeOnly[];
extern char txtAltML[];
extern char txtUseCoarseMesh[];
extern char txtOptimize[];
extern char txtPreTransf[];
extern char txtPreTransformation[];
extern char txtG1[];
extern char txtG2[];
extern char txtG1quasiG2[];
extern char txtRender[];
extern char txtCamera[];
extern char txtInterrupt[];
extern char txtLight[];
extern char txtC[];
extern char txtD_shape_func[];
extern char txtGaussian[];
extern char txtMean[];
extern char txtIsophotes[];
extern char txtReflection[];
extern char txtHighlight[];
extern char txtSections[];
extern char txtReflectionFrame[];
extern char txtHighlightFrame[];
extern char txtSectionsNormal[];
extern char txtShadows[];
extern char txtAntialias[];
extern char txtShapeFunction[];
extern char txtLight0[];
extern char txtLight1[];
extern char txtLight2[];
extern char txtLight3[];
extern char txtAmbient[];
extern char txtGeneral[];
extern char txtBlendingG1[];
extern char txtBlendingG2[];
extern char txtSwept[];
extern char txtSpherical[];
extern char txtLofted[];
extern char txtBiharmonic[];
extern char txtTriharmonic[];
extern char txtClamped[];
extern char txtColour[];
extern char txtSelectVertex[];
extern char txtSelectEdge[];
extern char txtOnLine[];
extern char txtOnPlane[];
extern char txtOnSphere[];
extern char txtOnCylinder[];
extern char txtRefPoint[];
extern char txtRefVector[];
extern char txtRadius[];
extern char txtAngleDeg[];
extern char txtX[];
extern char txtY[];
extern char txtZ[];
extern char txtProject[];
extern char txtLogIt[];
extern char txtNPThreads[];
extern char txtSpecials[];
extern char txtQuestionmark[];
extern char txtReset[];
extern char txtPosition[];
extern char txtF[];
extern char txtEquator[];
extern char txtMeridian[];
extern char txtSphericalProduct[];
extern char txtG1BlendingPatch[];
extern char txtG2BlendingPatch[];
extern char txtGetSubmesh[];
extern char txtPipeDiameter[];
extern char txtMengerCurv[];
extern char txtExponent[];
extern char txtPenaltyParam[];
extern char txtPenaltyOpt[];
extern char txtMoveMany[];

extern char *InfoMsg[];

extern char MsgReconsider[];

extern char ErrorMsgCannotOpen[];
extern char ErrorMsgCannotSave[];
extern char ErrorMsgIncorrectFilename[];
extern char ErrorMsgCannotRemoveVertex[];
extern char ErrorMsgCannotContractEdge[];
extern char ErrorMsgCannotGlueEdges[];
extern char ErrorMsgCannotGlueLoops[];
extern char ErrorMsgCannotSealHole[];
extern char ErrorMsgCannotRemoveFacet[];
extern char ErrorMsgCannotDoubleFacetEdges[];
extern char ErrorMsgMeshIntegrity[];
extern char ErrorMsgCannotClose[];
extern char ErrorMsgCannotRender[];
extern char ErrorMsgCannotLaunchAChild[];
extern char ErrorMsgChildProcessBusy[];
extern char ErrorMsgChildProcessFailure[];
extern char ErrorMsgCannotDoIt[];
extern char ErrorMsgOptimizationFailure[];
extern char ErrorMsgNotImplemented[];
extern char ErrorMsgNothingToMark[];
extern char ErrorMsgMeshCannotBeOptimized[];
extern char ErrorMessageCannotCangeDegree[];
extern char ErrorMsgCannotFlip[];
extern char ErrorMsgCannotFindObject[];
extern char ErrorMsgCouldNotExtractMesh[];
extern char ErrorMsgMustBeBoundaryEdge[];
extern char ErrorMsgNoValidPairOfVertices[];
extern char ErrorMsgCurveMustBeCubicAndClosed[];
extern char ErrorMsgCurveMustBeNonRationalAnd3D[];

