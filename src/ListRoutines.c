/* ------------------------------------------------------------------------------------
//  ListRoutines.c -- General list routines.  Adapted from an idea from Ben Halpern.
//  Written by JAT 1995;  
// ------------------------------------------------------------------------------------
*/
#include <stdlib.h>
#include <string.h>

#include "ListRoutines.h"

/*------------------------------------------------------------------------------------
//      CREATE NEW LIST -- Allocate and initialize a list objects of size sizeofobject.  
//                                         For tfoolist *foolist, use: foolist = CreateNewList(sizeof(tfoo), 10 , 10); 
//      A pointer to the list is returned.  If the list could not be created then NULL is returned.
*/
void *CreateNewList(INT_4 sizeofobject, INT_4 initialNum, INT_4 growNum) {

        tlist *list = nil;
        
        list = (tlist *) malloc(sizeof(tlist));
        if ( nil == list ) return nil;
        
        list->numObjects   = 0;
        list->limit        = initialNum;
        list->sizeofobject = sizeofobject;
        list->growNum      = growNum;
        list->object       = (char*) malloc( (size_t)(initialNum * sizeofobject) );
                
        if ( nil == list->object ) {
                free( list );
                return nil;
        }       

        return list;
}
/*------------------------------------------------------------------------------------
//  COPY LIST -- Make a copy of the list.
*/
void *CopyList(void *voidlist)
{
    tlist *inList  = nil;
        tlist *outList = nil;
    
    inList = (tlist *) voidlist;
        
        outList = (tlist *) CreateNewList(inList->sizeofobject, inList->numObjects, inList->growNum);
        if ( nil == outList ) return nil;

        outList->numObjects   = inList->numObjects;

        memcpy(&(outList->object[0]), &(inList->object[0]), 
                   (size_t)(outList->sizeofobject * outList->numObjects));
        
        return outList;
}
/*------------------------------------------------------------------------------------
//  ADD TO LIST -- Add an object of the appropriate size to the list. For tfoolist 
//                                 *foolist which was obtained with CreateNewList, and tfoo *foo, 
//  use: if ( !AddToList(foo, foolist) ) goto finishUp; The function returns true if
//  the object was successfully added to the list and false if it was not.
*/
char AddToList(void *voidobject, void *voidlist) {

        char     *newobjects;   
    toObject *object;
    tlist    *list;
    INT_4     sizeoflist;
    
    object  =   (toObject *) voidobject;
    list    =   (tlist *)   voidlist;
    sizeoflist   =   list->numObjects * list->sizeofobject;

    if (list->numObjects >= list->limit) {
                /* Boost the size of the list if we are going to overflow. */

        newobjects =  (char*) malloc((list->limit + list->growNum) * list->sizeofobject);
                if ( newobjects ) {
                memcpy(newobjects, list->object, (size_t)sizeoflist);
                
                list->limit += list->growNum;
                
                free( list->object );
                
                list->object = newobjects;
                }
                else return false;
    }
    
    memcpy(&(list->object[0]) + (size_t)sizeoflist, object, (size_t)(list->sizeofobject));
    
    list->numObjects++;
        return true;
}
/*------------------------------------------------------------------------------------
//  REMOVE FROM LIST -- Remove the object of the given index from the list.
*/
void RemoveFromList(INT_4 index, void *voidlist) {

        tlist *list;
        INT_4 sizeofobjectstomove;
        INT_4 sourceoffset, targetoffset;
        
    list    =   (tlist *)   voidlist;

        targetoffset = index * list->sizeofobject;
        sourceoffset = targetoffset + list->sizeofobject;
        sizeofobjectstomove = list->numObjects * list->sizeofobject - sourceoffset;
        
        memmove(list->object + targetoffset, list->object + sourceoffset, (size_t)sizeofobjectstomove);
        
        list->numObjects--;
}
/*------------------------------------------------------------------------------------
//  TRIM LIST -- Free up unused memory and reset limit accordingly.
*/
void TrimList(void *voidlist) {

        tlist *list;
        INT_4  sizeoflist;
        
        list = (tlist *) voidlist;
        
        sizeoflist = list->numObjects * list->sizeofobject;
        list->object = (char *) realloc(list->object, (size_t)sizeoflist );

        list->limit = list->numObjects;
}
/*------------------------------------------------------------------------------------
//  DISPOSE LIST -- Dispose of the list.
*/
void DisposeList(void *voiddeadlist) {

        tlist *deadlist;
        

        if ( nil == voiddeadlist ) return; /* Don't try to dispose of it if it doesn't exist.
                                                                                   Nasty things could happen. */
        deadlist = (tlist *) voiddeadlist;
        
        free( deadlist->object );
        free( deadlist );
        
}
