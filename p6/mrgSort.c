#include <stdlib.h>
#ifdef DEBUG
   #include <stdio.h>
   #include <assert.h>
#endif
#include "word.h"
// #define THUMB 1

#ifdef THUMB
void *mrgChains
(
   void *lb,
   void *rb,
   unsigned int gNxt,
   unsigned int gInt,
   void *endL,
   void *endR,
   void *ML
);

void *mrgSortR
(
   unsigned int N,
   void *bp,  
   unsigned long gNxt,
   unsigned long gI,
   void *ML, 
   void *UB 
);

#else
word_t *mrgChainsC
(
   word_t *lb,      /* address of left base node */
   word_t *rb,      /* address of right base node */
   word_t *endL,    //addr of end node of lb
   word_t *endR,    //addr of end node of rb
   word_t **ML      //OUT: points to end of merged list
)
/*
 * Given two sorted lists started by lb & rb.
 * NOTE: In my implementation, input and output lists NULL terminated.
 * RETURNS: base ptr for merged and sorted list, sets ML to last node
 */
{
   /*ptrs for base of new list, current node, and next node*/
   word_t *newListHead, *curr;
   int lbVal, rbVal;   //holds top val of lb and rb
   int endValL, endValR;   //holds last val of lb and rb
   //ASSUMES: lb,rb,endL,endR non-NULL due to mrgSortR
   lbVal = lb->len;
   rbVal = rb->len;
   endValL = endL->len;
   endValR = endR->len;
   //if greatest val in lb <= smallest val in rb, append rb to lb
   if(endValL <= rbVal)
   {
      *ML = endR;
      endL->next=rb;
      return lb;
   }
   //if greatest val in rb <= smallest val in lb, append lb to rb
   else if (endValR <= lbVal)
   {
      *ML = endL;
      endR->next=lb;
      return rb;
   } 
   /*initialize first item in new list by comparing left and right*/
   if(lbVal < rbVal)
   {
      newListHead = lb; /*initialize the base ptr for new list*/
      lb = lb->next;    /*increment the top of the left buffer*/
      if(!lb)           //if lb empty,
      {
         *ML = endR;
         newListHead->next=rb;   //append rb
         return newListHead;     //return newList
      }
      lbVal = lb->len;  //otherwise, update lbVal
   }
   else
   {
      newListHead = rb; /*initialize the base ptr for new list*/
      rb = rb->next;    /*increment the top of the right buffer*/
      if(!rb)           //if rb empty,
      {
         *ML = endL;
         newListHead->next=lb;   //append lb
         return newListHead;     //return newList
      }
      rbVal = rb->len;  //otherwise, update rbVal
   }
   curr = newListHead;
   /*comparison loop, using both lb and rb*/
   while(1)
   {
      /*compare the first item in the left to the first item in the right*/
      if(lbVal < rbVal)
      {
         curr->next = lb;     //next node is lb
         curr = lb;           //curr=curr->next
         lb = lb->next;       //increment lb
         if(!lb)              //if lb empty,
         {
            *ML = endR;
            curr->next=rb;       //append rb
            return newListHead;
         }
         lbVal = lb->len;  //otherwise, update lbVal
      }
      /*otherwise, set the first item in right as the first in new list*/
      else
      {
         curr->next = rb;     //next node is rb
         curr = rb;           //curr=curr->next
         rb = rb->next;       //increment rb
         if(!rb)              //if rb empty,
         {
            *ML = endL;
            curr->next=lb;       //append lb
            return newListHead;
         }
         rbVal = rb->len;  //otherwise, update rbVal
      }
   } //end while
}

/*
 * Recursive function:
 */
word_t *mrgSortRC   /* RETURNS: base of greatest-to-least ->len sorted list */
(                   /*          using merge sort, so O( N log2(N) )         */
   unsigned int N,  /* number of nodes in list ub */
   word_t *bp,      /* base ptr for N nodes, not necessarily NULL-terminated */
   unsigned long gNxt,
   unsigned long gI,
   word_t **ML,     /* OUT: I used this to point last node in merged list */
   word_t **UB      /* OUT: *UB set to address of n+1 node, which is unsorted */
)
/*
 * On entry ub is an unsorted list, which will be destroyed to produce sort.
 * NOTE: Input lists may not be NULL terminated, but output lists will be.
 * RETURNS: base ptr of list sorted by ->len.
 */
{
   if(N == 1) //base case, sorted
   {
      /*increment *UB so it is always pointing to the next unsorted pos*/
      *UB = bp->next;
      //last node in list is only node in list, bp
      *ML = bp;
      //null terminate
      bp->next = NULL;
      return bp;
   }
   //bp is always non-NULL and N>0 thanks to wrapper function
   word_t *lb, *rb;
   word_t *endL, *endR;  //ptrs to last nodes of lb and rb
   /*initialize lengths of left and right buffers for splitting*/
   unsigned int leftHalf=(N>>1), rightHalf = N-leftHalf;
   /*mrgSort on left nodes*/
   lb = mrgSortRC(leftHalf,bp,gNxt,gI,&endL,UB);
   //rb set to next unsorted spot
   //endL set to last node in lb
   /*mrgSort on right nodes, utilizing the new rb value from prev mrgSort*/
   rb = mrgSortRC(rightHalf,*UB,gNxt,gI,&endR,UB);
   //UB updated to next unsorted spot
   //endR set to last node in rb
   
   //return mergeLists on left and right
//   return mrgChains(lb,rb,gNxt,gI,endL,endR,ML);
   return mrgChainsC(lb,rb,endL,endR,ML);
}
#endif //ifdef THUMB

/*
 * Wrapper function
 */
word_t *mrgSort(word_t *bp, unsigned int gNxt, unsigned int gI)
/*
 * NOTE: mrgSortFast is not recursive: the optimized mrgSort above handles
 *       the recursion while passing extra information (not of interest to
 *       top level) in order to avoid redundant memory access.
 * RETURNS: address of base node of list sorted least-to-greatest on ->len.
 *          ub's nodes are used to create list, so unsorted list will be
 *          destroyed by sort.
 */
{
   if(bp) /*ensure there's at least 1 element*/
   {
      register int N=1, prevLen=bp->len, sorted=1;
      word_t *curr, *UB=NULL, *ML=NULL;

      /*find length of list and determine whether it's sorted already*/
      for(curr=bp->next; curr; curr=curr->next)
      {
         register int currLen = curr->len;
         sorted &= (prevLen <= currLen); /*remains 1 if sorted, 0 if not*/
         /*update prevLen before incrementing curr*/
         prevLen = currLen;
         N++;
      }
      /*if found to be sorted, return*/
      if(sorted)
         return(bp);
      /*otherwise, sort*/
      #ifdef THUMB
      return mrgSortR(N,bp,gNxt,gI,&ML,&UB);
      #else
      return mrgSortRC(N, bp, gNxt, gI, &ML, &UB);
      #endif
   }
   return bp;
}
