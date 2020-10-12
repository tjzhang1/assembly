#include <stdlib.h>
#ifdef DEBUG
   #include <stdio.h>
   #include <assert.h>
#endif
#include "word.h"

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
   word_t *current_lb=lb, *current_rb=rb;
   /*ptrs for base of new list and the current end*/
   word_t *newListHead, *curr;
   /*used to count how many elements are left in lb and rb*/
   register unsigned int lbLen=nL, rbLen=nR;

   if(!rb)
      return lb;
   if(!lb)
      return rb;

   /*initialize first item in new list by comparing left and right*/
   if(lb->len <= rb->len)
   {
      newListHead = lb; /*initialize the base ptr for new list*/
      current_lb = lb->next; /*increment the top of the left buffer*/
      lbLen--; /*reduce the num of elements left in left buffer by 1*/
   }
   else
   {
      newListHead = rb; /*initialize the base ptr for new list*/
      current_rb = rb->next; /*increment the top of the right buffer*/
      rbLen--; /*reduce the num of elements left in right buffer by 1*/
   }
   /*comparison loop, using both lb and rb*/
   for(curr=newListHead; lbLen&&rbLen; curr=curr->next)
   {
      /*compare the first item in the left to the first item in the right*/
      if(current_lb->len <= current_rb->len)
      {
         curr->next = current_lb;
         current_lb = current_lb->next;
         lbLen--;
      }
      /*otherwise, set the first item in right as the first in new list*/
      else
      {
         curr->next = current_rb;
         current_rb = current_rb->next;
         rbLen--;
      }
   }
   /*if the left still has stuff, add the rest to end of curr*/
   while( lbLen )
   {
      curr->next = current_lb;
      current_lb = current_lb->next;
      curr=curr->next;
      lbLen--;
   }
   /*if the right still has stuff, add the rest to end of curr*/
   while( rbLen )
   {
      curr->next = current_rb;
      current_rb = current_rb->next;
      curr=curr->next;
      rbLen--;
   }
   /*all lists are empty, NULL terminate the list*/
   curr->next=NULL;

   return newListHead;
}

word_t *mrgSort     /* RETURNS: base of greatest-to-least ->len sorted list */
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
      /*return ub because a 1 length list is sorted*/
      return ub;
   }

   word_t *lb, *rb;
   /*initialize lengths of left and right buffers for splitting*/
   register int leftHalf=(N>>1), rightHalf = N-leftHalf;
   /*mrgSort on left nodes*/
   lb = mrgSort(leftHalf, ub, UB); /*this will increment *UB ptr*/
   /*mrgSort on right nodes, utilizing the new UB value from prev mrgSort*/
   rb = mrgSort(rightHalf, *UB, UB);
   /*return mergeLists on left and right*/
   return mergeLists(leftHalf, rightHalf, lb, rb);

}

word_t *mrgSortFast(word_t *ub) /* required ABI, wraps recursive mrgSort */
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
      register int N=1, prevLen=ub->len, unsorted=0;
      word_t *curr, *UB=ub;

      /*find length of list and determine whether it's sorted already*/
      for(curr=ub->next; curr; curr=curr->next)
      {
         register int currLen = curr->len;
         unsorted += (prevLen > currLen); /*remains 0 if sorted, > 0 if not*/
         /*update prevLen before incrementing curr*/
         prevLen = currLen;
         N++;
      }
      /*if found to be sorted, return*/
      if(!unsorted)
         return(ub);
      /*otherwise, sort*/
      return mrgSort(N, ub, &UB);
   }
   return ub;
}
