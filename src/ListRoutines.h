/* ------------------------------------------------------------------------------------
//  ListRoutines.h -- definitions and prototypes for ListRoutines.c
//               
// ------------------------------------------------------------------------------------
*/
#pragma once

#ifndef __LISTROUTINES__
#define __LISTROUTINES__

#include "LutefiskDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#ifndef nil
#define nil 0
#endif

/* Generic list structues */

typedef char toObject;

typedef struct {
        INT_4     numObjects;
        INT_4     limit;
        INT_4     sizeofobject;
        INT_4     growNum;
        toObject *object;
} tlist;

typedef struct {
        INT_4  numObjects;
        INT_4  limit;
        INT_4  sizeofobject;
        INT_4  growNum;
        INT_4 *entry;
} tINT_2list;

/* Public function prototypes ---~---~---~---~---~---~---~---~---~---~---~---~---~ */
extern void *CreateNewList(INT_4 sizeofobject, INT_4 initialNum, INT_4 growNum);
extern void *CopyList(void *voidlist);
extern char AddToList(void *voidobject, void *voidlist);
extern void RemoveFromList(INT_4 index, void *voidlist);
extern void TrimList(void *voidlist);
extern void DisposeList(void *voiddeadlist);


#ifdef __cplusplus
}
#endif

#endif /* __LISTROUTINES__ */   
