
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_MEM_BLOCKS 32

void MemBlockListInit ( void );
void MemBlockRegister ( void *ptr, boolean alloc );
void MemBlockFreeAll ( void );

