#include <stdlib.h>
#ifdef DEBUG
   #include <stdio.h>
   #include <assert.h>
#endif
#include "word.h"

#ifdef THUMB_CHAINS
void *mrgChains
(
   void *lb,
   void *rb,
   unsigned int gNxt,
   unsigned int gInt,
   void *ML,
   void *endL,
   void *endR
);
#else
word_t *mrgChainsC
(
   word_t *lb,      // address of left base node 
   word_t *rb,      // address of right base node
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
   //ASSUME: lb will always run out before rb, saving us extra checks
   //   This always occors if max lb <= max rb
   //if max lb > max rb, swap lb and rb to make this assumption true
   if (endValL > endValR)
   {
      word_t *tmp;
      //last node of merged list is endL now
      *ML = endL;
      //swap lb and rb
      tmp = lb;
      lb = rb;
      rb = tmp;
      //also swap lbLen and rbLen
      lbVal ^= rbVal;
      rbVal ^= lbVal;
      lbVal ^= rbVal;
   }
   else
      //the last node of merged list is endR
      *ML = endR; 
   /*initialize first item in new list by comparing left and right*/
   if(lbVal <= rbVal)
   {
      newListHead = lb; /*initialize the base ptr for new list*/
      lb = lb->next;    /*increment the top of the left buffer*/
      if(!lb)           //if lb empty,
      {
         newListHead->next=rb;   //append rb
         return newListHead;     //return newList
      }
      lbVal = lb->len;  //otherwise, update lbVal
   }
   else
   {
      newListHead = rb; //initialize the base ptr for new list
      rb = rb->next;    //increment the top of the right buffer
      //don't need to check if rb->next is empty, lb empty before rb
      rbVal = rb->len;  //update rbVal
   }
   curr = newListHead;
   /*comparison loop, using both lb and rb*/
   while(1)
   {
      //compare the first item in the left to the first item in the right
      if(lbVal <= rbVal)
      {
         curr->next = lb;     //next node is lb
         curr = lb;           //curr=curr->next
         lb = lb->next;       //increment lb
         if(!lb)              //if lb empty,
         {
            curr->next=rb;       //append rb
            return newListHead;
         }
         lbVal = lb->len;  //otherwise, update lbVal
         //continue left: don't need to set the next pointer b/c curr->next
         //already points to the next item in lb
         while(lbVal <= rbVal)
         {
            curr = lb;
            lb = lb->next;
            if(!lb)
            {
               curr->next=rb;  //append rb
               return newListHead;
            }
            lbVal = lb->len;
         }
      }
      //otherwise, set the first item in right as the next
      else
      {
         curr->next = rb;     //next node is rb
         curr = rb;           //curr=curr->next
         rb = rb->next;       //increment rb
         //don't need to check if rb->next is empty, lb empty before rb
         rbVal = rb->len;     // update rbVal
         //continue right: don't need to set the next pointer b/c curr-next
         //already points to the next item in rb
         //can't be <= or else there's a chance rb runs out of nodes first
         while(rbVal < lbVal)  
         {
            curr = rb;
            rb = rb->next;
            //don't need to check if rb->next is empty, lb empty before rb
            rbVal = rb->len;
         }
      }
   } //end while
}
#endif

/*
 * Recursive function:
 */
#ifdef THUMB_SORTR
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
   else if (N == 2)  //N = 2 case
   {
      word_t *n0=bp, *n1=bp->next; //node 1 and node 2
      int n0Val=n0->len, n1Val=n1->len;

      *UB = n1->next;
      if(n0Val <= n1Val)  //already sorted
      {
         *ML = n1;
         n1->next = NULL;
         return n0;
      }
      //if not sorted, then swap n0 and n1
      n1->next = n0;
      n0->next = NULL;
      *ML = n0;
      return n1;
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
   #ifdef THUMB_CHAINS
   return mrgChains(lb,rb,gNxt,gI,ML,endL,endR);
   #else
   return mrgChainsC(lb,rb,endL,endR,ML);
   #endif
}
#endif 

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
      #ifdef THUMB_SORTR
      return mrgSortR(N,bp,gNxt,gI,&ML,&UB);
      #else
      return mrgSortRC(N, bp, gNxt, gI, &ML, &UB);
      #endif
   }
   return bp;
}
