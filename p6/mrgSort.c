#include <stdlib.h>
#ifdef DEBUG
   #include <stdio.h>
   #include <assert.h>
#endif
#include "word.h"
#define CHAINS 1

void *mrgChains4(void *lb, void *rb, unsigned int gNxt, unsigned int gInt);

word_t *mergeLists
(
   unsigned int nL, /* # of nodes in unsorted list lb (left base) */
   unsigned int nR, /* # of nodes in unsorted list rb (Right base) */
   word_t *lb,      /* address of left base node */
   word_t *rb       /* address of right base node */
)
/*
 * Given two sorted lists started by lb & rb.
 * NOTE: In my implementation, input and output lists NULL terminated.
 * RETURNS: base ptr for merged  and sorted list
 */
{
   /*ptrs to top of left and right buffers*/
   /*ptrs for base of new list, current node, and next node*/
   word_t *newListHead, *curr;
   /*used to count how many elements are left in lb and rb*/
   register unsigned int lbLen=nL, rbLen=nR;
   register unsigned int lbVal, rbVal;

   if(!rb)
      return lb;
   if(!lb)
      return rb;
   lbVal = lb->len;
   rbVal = rb->len;
   /*initialize first item in new list by comparing left and right*/
   if(lbVal < rbVal)
   {
      newListHead = lb; /*initialize the base ptr for new list*/
      lb = lb->next;    /*increment the top of the left buffer*/
      lbLen--;          /*reduce the num of elements left in left buffer by 1*/
      lbVal = (lbLen) ? lb->len:lbVal;
   }
   else
   {
      newListHead = rb; /*initialize the base ptr for new list*/
      rb = rb->next; /*increment the top of the right buffer*/
      rbLen--; /*reduce the num of elements left in right buffer by 1*/
      rbVal = (rbLen) ? rb->len:rbVal;
   }
   /*comparison loop, using both lb and rb*/
   for(curr=newListHead; lbLen&&rbLen; curr=curr->next)
   {
      /*compare the first item in the left to the first item in the right*/
      if(lbVal < rbVal)
      {
         curr->next = lb;
         lb = lb->next;
         lbLen--;
         lbVal = (lbLen) ? lb->len:lbVal;
      }
      /*otherwise, set the first item in right as the first in new list*/
      else
      {
         curr->next = rb;
         rb = rb->next;
         rbLen--;
         rbVal = (rbLen) ? rb->len:rbVal;
      }
   }
   //if the left still has stuff, add the rest to end of curr

   while( lbLen )
   {
      curr->next=lb;
      curr=lb;
      lb = lb->next;
      lbLen--;
   }
   //if the right still has stuff, add the rest to end of curr
   while( rbLen )
   {
      curr->next=rb;
      curr=rb;
      rb = rb->next;
      rbLen--;
   }
   //all lists are empty, NULL terminate the list
   curr->next=NULL;

   return newListHead;
}

word_t *mrgChains
(
   unsigned int a,
   unsigned int b,
   word_t *lb,
   word_t *rb
)
/*
 * Given two sorted lists started by lb & rb.
 * NOTE: In my implementation, input and output lists NULL terminated.
 * RETURNS: base ptr for merged  and sorted list
 */
{
   /*ptrs to top of left and right buffers*/
   /*ptrs for base of new list, current node, and next node*/
   word_t *newListHead, *curr;
   /*used to count how many elements are left in lb and rb*/
   register unsigned int lbVal, rbVal;

   if(!rb)
      return lb;
   if(!lb)
      return rb;
   lbVal = lb->len;
   rbVal = rb->len;
   /*initialize first item in new list by comparing left and right*/
   if(lbVal < rbVal)
   {
      newListHead = lb; /*initialize the base ptr for new list*/
      lb = lb->next;    /*increment the top of the left buffer*/
      lbVal = (lb) ? lb->len:lbVal;
   }
   else
   {
      newListHead = rb; /*initialize the base ptr for new list*/
      rb = rb->next; /*increment the top of the right buffer*/
      rbVal = (rb) ? rb->len:rbVal;
   }
   /*comparison loop, using both lb and rb*/
   for(curr=newListHead; lb&&rb; curr=curr->next)
   {
      /*compare the first item in the left to the first item in the right*/
      if(lbVal < rbVal)
      {
         curr->next = lb;
         lb = lb->next;
         lbVal = (lb) ? lb->len:lbVal;
      }
      /*otherwise, set the first item in right as the first in new list*/
      else
      {
         curr->next = rb;
         rb = rb->next;
         rbVal = (rb) ? rb->len:rbVal;
      }
   }

   if(lb)
      curr->next=lb;
   else
      curr->next=rb;

   return newListHead;
}

word_t *mrgSortR    /* RETURNS: base of greatest-to-least ->len sorted list */
(                   /*          using merge sort, so O( N log2(N) )         */
   unsigned int N,  /* number of nodes in list ub */
   word_t *ub,      /* base ptr for N nodes, not necessarily NULL-terminated */
   word_t **UB      /* OUT: *UB set to address of n+1 node, which is unsorted */
)
/*
 * On entry ub is an unsorted list, which will be destroyed to produce sort.
 * NOTE: Input lists may not be NULL terminated, but output lists will be.
 * RETURNS: base ptr of list sorted by ->len.
 */
{
   /*
   Pseudocode
   1. check if N is 1; if so increment *UB=ub->next and return ub
   2. otherwise call mrgSort on N/2 nodes, save as left
   3. call mrgSort on N-N/2 nodes, save as right
   4. return mergeLists on left and right
   */
   if(N == 1) /*base case*/
   {
      /*increment *UB so it is always pointing to the next unsorted pos*/
      *UB = ub->next;
      #ifdef CHAINS
         ub->next = NULL;
      #endif
      /*return ub because a 1 length list is sorted*/
      return ub;
   }
   //ub is always non-NULL and N>0 thanks to wrapper function
   word_t *lb, *rb;
   /*initialize lengths of left and right buffers for splitting*/
   register int leftHalf=(N>>1), rightHalf = N-leftHalf;
   /*mrgSort on left nodes*/
   lb = mrgSortR(leftHalf, ub, &rb); /*this will set rb to next unsorted spot*/
   /*mrgSort on right nodes, utilizing the new UB value from prev mrgSort*/
   rb = mrgSortR(rightHalf, rb, UB);
   /*return mergeLists on left and right*/
   #if CHAINS
//      return mrgChains(leftHalf, rightHalf, lb, rb);
      return mrgChains4(lb,rb,0,0);
   #else
      return mergeLists(leftHalf,rightHalf,lb,rb);
   #endif
}

word_t *mrgSort(word_t *ub, int a, int b) /* required ABI, wraps recursive mrgSort */
/*
 * NOTE: mrgSortFast is not recursive: the optimized mrgSort above handles
 *       the recursion while passing extra information (not of interest to
 *       top level) in order to avoid redundant memory access.
 * RETURNS: address of base node of list sorted least-to-greatest on ->len.
 *          ub's nodes are used to create list, so unsorted list will be
 *          destroyed by sort.
 */
{
   /*pseudocode
   1. Traverse list to find length for N
      - determine if sorted already along the way
   2. Assign *UB to ub
   3. mrgSort(N, ub, UB)
   */
   if(ub) /*ensure there's at least 1 element*/
   {
      register int N=1, prevLen=ub->len, sorted=1;
      word_t *curr, *UB=ub;

      /*find length of list and determine whether it's sorted already*/
      for(curr=ub->next; curr; curr=curr->next)
      {
         register int currLen = curr->len;
         sorted &= (prevLen <= currLen); /*remains 1 if sorted, 0 if not*/
         /*update prevLen before incrementing curr*/
         prevLen = currLen;
         N++;
      }
      /*if found to be sorted, return*/
      if(sorted)
         return(ub);
      /*otherwise, sort*/
      return mrgSortR(N, ub, &UB);
   }
   return ub;
}
