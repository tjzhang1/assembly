#include <stdlib.h>
#include <assert.h>
#include "word.h"

word_t *mergeLists(word_t *lb, word_t *rb)
/*
 * Given two (left & right) sorted lists started by lb & rb
 * RETURNS: base ptr for merged  and sorted list
 * NOTE: both lb & rb are lost, since their nodes are merged into return list
 */
{
   word_t *newListHead, *current_nl;
   word_t *current_lb, *current_rb;
   /*error checking if one or both lb and rb are NULL*/
   if(lb && !rb)
      return lb;
   else if(!lb && rb)
      return rb;
   else if(!lb && !rb)
      return NULL;
   /*initialize variables to hold current positions in lb and rb*/
   current_lb = lb;
   current_rb = rb;
   /*initialize newListHead*/
   /*compare the first item in the left to the first item in the right*/
   if(current_lb->len <= current_rb->len)
   {
      newListHead = current_lb;
      current_lb = current_lb->next;
   }
   else /*otherwise, set the first item in right as the first in new list*/
   {
      newListHead = current_rb;
      current_rb = current_rb->next;
   }
   current_nl = newListHead;

   /*loop until one or both lists are empty*/
   while(current_lb && current_rb)
   {
      /*compare the first item in the left to the first item in the right*/
      if(current_lb->len <= current_rb->len)
      {
         current_nl->next = current_lb;
         current_lb = current_lb->next;
      }
      else /*otherwise, set the first item in right as the first in new list*/
      {
         current_nl->next = current_rb;
         current_rb = current_rb->next;
      }
      /*increment current_nl pointer*/
      current_nl = current_nl->next;
   }

   /*if the left still has stuff, add the rest to end of current_nl*/
   if(current_lb && !current_rb)
      current_nl->next = current_lb;
   /*if the right still has stuff, add the rest to end of current_nl*/
   else if(!current_lb && current_rb)
      current_nl->next = current_rb;
   /*otherwise, left and right are both empty, end list with NULL*/
   else
      current_nl->next = NULL;

   return(newListHead);
}

word_t *mrgSortSlow /*RET: base ptr of list sorted least-to-greatest by ->len.*/
(
   word_t *ub /* base of unsorted list, destroyed by sort */
)
/*
 * On entry ub is an unsorted list, which will be destroyed to produce sort.
 * RETURNS: base ptr of list sorted least-to-greatest by ->len.
 * NOTE: for this API must be used for all recursion, and no external variables
 *       can be used (mrgSortFast will allow extending recursive params).
 */
{
   /*
   Pseudocode
   1. find length of linked list
     - can also check if sorted
     - can simultaneously find center
   2. if linked list is sorted, return it
   3. split left and right half
   4. mergeSort on left
   5. mergeSort on right
   6. return mergeLists on left and right half
   */
   if(ub)
   {
      if(ub->next)
      {
         register int len, unsorted, prevLen;
         word_t *left, *right; /*pointers to hold left and right halves*/
         word_t *curr; /*pointers to traverse linked list*/
         /*pointers to hold the center and the node before*/
         word_t *center, *center_prev;

         /*assume sorted until proven wrong*/
         unsorted = 0;
         /*initialize center pointers*/
         center_prev = NULL;
         center = ub;
         prevLen = ub->len;
         len = 1;
         for (curr=ub->next; curr; curr = curr->next)
         {
            register int currLen = curr->len;
            /*if previous is <= current (sorted), the value of unsorted remains 0*/
            /*otherwise, the value of unsorted increments*/
            unsorted += (prevLen > currLen);
            prevLen = curr->len;
            len++;
            /*increment center once every other time (i.e. len is even);
             len is even when the least significant digit is 0*/
            if(!(len & 0x01))
            {
               center_prev = center;
               center = center->next;
            }
         }
         /*if linked list is sorted, return ub*/
         if(!unsorted)
            return ub;
         /*split left and right halves*/
         left = ub;
         right = center;
         center_prev->next = NULL; /*close the left half with NULL*/

         /*merge sort on left and right*/
         left = mrgSortSlow(left);
         right = mrgSortSlow(right);

         /*merge lists*/
         return mergeLists(left, right);
      }
   }
   return ub;
}
