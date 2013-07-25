
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

typedef struct kl_link {
    struct kl_joint *firstjoint;
    trans3d M;
    int    id;
    void (*drawlink) ( struct kl_link *lnk, trans3d *tr );
  } kl_link;

typedef struct kl_joint {
    struct kl_joint *next;
    struct kl_link *lnk, *up;
    trans3d L, R;
    int    id;
    void (*articulate) ( struct kl_joint *jnt );
  } kl_joint;


kl_link *NewLink ( int id, void (*drawlink) ( struct kl_link *lnk, trans3d *tr ) );
kl_joint *NewJoint ( int id, kl_link *first, kl_link *second,
                     void (*articulate) ( struct kl_joint *jnt ) );

void DisplayLinkage ( kl_joint *jnt, trans3d *m );

