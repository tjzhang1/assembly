#include <stdlib.h>
#include <assert.h>

#define uint unsigned int
/*
 * These macros completely self-explanatory, and clear by noob inspection.
 * In case not, they allow you to access (ld:load, st:store)  either pointers
 * or integers stored at any effective address (EA=ptr+offset) in memory.
 */
#define ldPtr(p_, off_) ( *((int**)(((char*)(p_))+(off_))) )
#define stPtr(p_, off_, v_) *((int**)(((char*)(p_))+(off_))) = (int*)(v_);
#define ldInt(p_, off_) ( *((int*)(((char*)(p_))+(off_))) )

void *mrgChains(void *lb, void *rb, uint gNxt, uint gI); // merge sorted lists

void *mrgSortR/* RETURNS: possibly changed base ptr of sorted list */
(
   uint N,    /* Number of nodes of list that starts at address in lb */
   void *lb,  /* on ENTRY unsorted list baseptr (later left side after split) */
   uint gNxt, /* gap in bytes to ->next ptr */
   uint gI,   /* gap in bytes to integer we are sorting on */
   void *UB   /* on EXIT: beginning of unsorted list */
)
/*
 * Recursively sorts N-length list lb, without assuming lb NULL terminated.
 * This is a generic function that works without knowing anything about
 * the structures except that they are singly-linked lists, with the next
 * ptr gNxt bytes from node beginning, and the integer sort value at offset gI.
 *
 * NOTE: it calls generic function mrgChains, which can merge the sorted lists
 *       in O(N) time,  N=length(left list)+length(right list).
 *
 * RETURNS: possibly changed base ptr of sorted list
 *    UB returns the base ptr of the unsorted nodes.
 */
{
   if (N == 1)  // Basis case makes lists NULL term, as expected by mrgChains
   {
      int **pp = UB;         // need to set UB, but can't deref void*
      *pp = ldPtr(lb, gNxt); // UB = lb->next */
      stPtr(lb, gNxt, NULL); // lb->next=NULL : lb 1-node (sorted) NULL-term 
      return(lb);            // return sorted in lb, unsorted *UB 
   }
   else if (N > 1)
   {
      unsigned int nL=N>>1, nR=N-nL;
      void *rb;
//
//    First call's unsorted base return tells me where my right search starts
//    but 2nd one says where me and all my children stopped searching, so
//    is passed back up to caller in shared space (local var in wrapper func!)
//
      lb = mrgSortR(nL, lb, gNxt, gI, &rb);// get nL-len sorted list wt base lb
      rb = mrgSortR(nR, rb, gNxt, gI, UB); // get nR-len sorted list wt base rb
      lb = mrgChains(lb, rb, gNxt, gI);    // merge sorted lists
      return(lb);                          // ret merged list, *UB set by rec
   }
   return(NULL);  // do not set *UB for N <= 0
}

void *mrgSort  // RETURNS: possibly changed baseptr for sorted N-len list
(
   unsigned int N,    // length of input list, ASSUMED correct
   void *bp,          // address of first node in list
   unsigned int gNxt, // gap in bytes from start of struct to next ptr
   unsigned int gI    // gap in bytes from start of struct to int to sort on
)
/*
 * This non-recursive wrapper allocates space for the unsorted base ptr (ub),
 * and then starts the recursion.
 */
{
   void *ub;
   return(mrgSortR(N, bp, gNxt, gI, &ub));
}
