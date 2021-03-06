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
   int nL,          //num nodes remaining in lb
   word_t *lb,      //address of left base node
   word_t *rb,      //address of right base node
   int lbVal,       //val of first node of lb
   int rbVal        //val of first node of rb
)
/*
 * Given two sorted lists started by lb & rb, returns merged and sorted list.
 * ASSUME: input and output lists NULL terminated
 * ASSUME: lb will run out of nodes before rb will
 * RETURNS: base ptr for merged and sorted list
 */
{
   /*ptrs for base of new list, current node, and next node*/
   word_t *newListHead, *curr;
   //Ptrs to preload lb->next and rb->next
   word_t *nextL=lb->next, *nextR=rb->next;
   /*initialize first item in new list by comparing left and right*/
   if(lbVal <= rbVal)
   {
      nL--;
      newListHead = lb; //initialize the base ptr for new list
      lb = nextL;       //increment the top of the left buffer
      lbVal = nextL->len;   //update lbVal
      //Don't need to check if lb empty because that would mean lb has only one
      //   element < minR. If lb has only one element, then minL=maxL < minR,
      //   which was already handled in mrgSortRC.
      nextL = nextL->next;  //begin loading next node into nextL
 //     nextLV = nextL->len;
   }
   else
   {
      newListHead = rb; //initialize the base ptr for new list
      rb = nextR;       //increment the top of the right buffer
      //don't need to check if rb is empty, lb empty before rb
      rbVal = nextR->len;  //update rbVal
      nextR = nextR->next;  //begin loading next node into nextR
   }
   curr = newListHead;
   //comparison loop, using both lb and rb
   while(1)
   {
      //compare the first item in the left to the first item in the right
      if(lbVal <= rbVal)
      {
         nL--;
         curr->next = lb;     //next node is lb
         curr = lb;           //curr=curr->next
         if(!nL)              //if lb empty,
         {
            curr->next=rb;       //append rb
            return newListHead;
         }
         lb = nextL;          //increment lb
         lbVal = nextL->len;  //update lbVal
         nextL = nextL->next; //begin loading next value of nextL
         //continue left: don't need to set the next pointer b/c curr->next
         //already points to the next item in lb
         while(lbVal <= rbVal)
         {
            nL--;
            //don't need to write curr->next
            curr = lb;
            if(!nL)
            {
               curr->next=rb;  //append rb
               return newListHead;
            }
            lb = nextL;
            lbVal = nextL->len;
            nextL = nextL->next;
         }
      }
      //otherwise, set the first item in right as the next
      else
      {
         curr->next = rb;     //next node is rb
         curr = rb;           //curr=curr->next
         rb = nextR;          //increment rb
         //don't need to check if rb is empty, lb empty before rb
         rbVal = nextR->len;  //update rbVal
         nextR = nextR->next; //increment nextR
         //continue right: don't need to set the next pointer b/c curr-next
         //   already points to the next item in rb
         //can't be <= or else there's a chance rb runs out of nodes first
         while(rbVal < lbVal)  
         {
            curr = rb;
            rb = nextR;
            //don't need to check if rb is empty, lb empty before rb
            rbVal = nextR->len;
            nextR = nextR->next;
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
   int minL, minR;       //holds top (min) val of lb and rb
   int maxL, maxR;       //holds last (max) val of lb and rb
   //initialize lengths of left and right buffers for splitting
   unsigned int leftHalf=N>>1, rightHalf = N-leftHalf;
   /*mrgSort on left nodes*/
   lb = mrgSortRC(leftHalf,bp,gNxt,gI,&endL,UB);
   //rb set to next unsorted spot
   //endL set to last node in lb
   /*mrgSort on right nodes, utilizing the new rb value from prev mrgSort*/
   rb = mrgSortRC(rightHalf,*UB,gNxt,gI,&endR,UB);
   //UB updated to next unsorted spot
   //endR set to last node in rb
   maxL = endL->len;
   minR = rb->len;
   minL = lb->len;
   maxR = endR->len;
   //if max val in lb <= min val in rb, append rb to lb
   if(maxL <= minR)
   {
      *ML = endR;
      endL->next=rb;
      return lb;
   }
   //if max val in rb <= smallest val in lb, append lb to rb
   else if (maxR <= minL)
   {
      *ML = endL;
      endR->next=lb;
      return rb;
   }
   //Make it so lb will always run out of nodes before rb. This always occurs
   //  if max lb <= max rb. Benefits: simplifying mrgChains loops by reducing
   //  num of arguments needed to pass in. Since lb runs out of
   //  nodes before rb, then only need to check how many nodes left in lb.
   //if max lb > max rb, swap lb and rb to make this assumption true
   if (maxL > maxR)
   {
      //last node of merged list is endL now
      *ML = endL;
      //call mrgChains with lb and rb swapped
      return mrgChainsC(rightHalf,rb,lb,minR,minL);
   }
   else
   {
      //the last node of merged list is endR
      *ML = endR;
      //call mrgChains with lb and rb in standard order
      return mrgChainsC(leftHalf,lb,rb,minL,minR);
   }
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
