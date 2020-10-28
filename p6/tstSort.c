#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

unsigned long myGap(void *p0, void *p1)
{
   size_t a0=(size_t)p0, a1=(size_t)p1;
   return((a0 <= a1) ? a1-a0 : a0-a1);
}

#ifndef WORD_ONLY
typedef struct items item_t;
struct items
{
   char name[32];       /* name of item */
   item_t *next;        /* ptr to next item in queue */
   float cost, weight;  /* mul by num for weight of all items */
   int i0, i1;          /* usually low & hi damage, but varies by item */
   unsigned int num;    /* number of this item being carried */
   unsigned int dsc;    /* if non-zero, extended description */
   unsigned int flg;    /* bitvec used to store info about item */
};

unsigned int tstItems(item_t *bp, int which)
{
   if (bp)
   {
      unsigned int nerr=0, i=1;
      int lastV = which ? bp->i1 : bp->i0;
      bp = bp->next;
      while(bp)
      {
         const int V = which ? bp->i1 : bp->i0;
         if (V < lastV)
         {
            fprintf(stderr, "ERROR on ->i%u, loc %u: val=%d, lastV=%d!\n",
                    which, i, V, lastV);
            nerr++;
         }
         i++;
         bp = bp->next;
      }
      return(nerr);
   }
   return(0);
}
#endif

#include "word.h"

int main(int nargs, char **args)
{
   unsigned int i, n;
   unsigned int gNxt, gL;
#ifndef WORD_ONLY
   unsigned int gN, g0, g1;
   item_t *ip, *ib=NULL;
#endif
   word_t *wp, *wb=NULL;
   unsigned long addrGap(void *a0, void *a1);
   void *mrgSort(void *bp, unsigned long gNxt, unsigned long gI);

   gNxt = myGap(wb, &(wb->next));
   gL = myGap(wb, &(wb->len));
   #ifndef WORD_ONLY
      gN = myGap(ib, &(ib->next));
      g0 = myGap(ib, &(ib->i0));
      g1 = myGap(ib, &(ib->i1));
   #endif
   for (n=1; n < nargs; n++)
   {
      unsigned int N, nerr=0;
      N = atoi(args[n]);
      fprintf(stderr, "TESTING N=%u case, word_t\n", N);
/*
 *    Allocate N-length lists of words and items, with random int vals
 *    make repeatable by length by seeding with length.
 */
      srand(N);
      for (i=0; i < N; i++)
      {
         wp = malloc(sizeof(word_t));
         assert(wp);
         wp->next = wb;
         wb = wp;
         #ifndef WORD_ONLY
            ip = malloc(sizeof(item_t));
            assert(ip);
            ip->next = ib;
            ib = ip;
         #endif
         #if 0
            wp->len = ip->i0 = ip->i1 = N-i;
         #else
            wp->len = rand();
            #ifndef WORD_ONLY
               ip->i0 = rand();
               ip->i1 = rand();
            #endif
         #endif
      }
      wp = wb;
      wb = mrgSort(wb, gNxt, gL);
      if (N > 1)
      {
         word_t *prevW=wb;
         int k;
         for (k=1,wp=wb->next; wp; wp = wp->next, k++)
         {
            if (prevW->len > wp->len)
            {
               fprintf(stderr, "node %u out-of-order, value=%d, prev=%d\n", 
                      k, wp->len, prevW->len);
               nerr++;
            }
            prevW = wp;
         }
         if (N != k)
         {
            if (N > k)
               fprintf(stderr, "SORT LOST %u NODES FROM LIST!\n\n", N-k);
            else if (k < N)
               fprintf(stderr, "SORT ADDED %u NODES!\n\n", k-N);
            nerr++;
         }
      }
      else if (N == 1)
          assert(wb==wp && !wb->next);
      else /* N == 0 */
         assert(wb == NULL);
      if (nerr)
      {
         fprintf(stderr, "Stopping testing, bad case N=%u words, NERR=%u\n",
                 N, nerr);
         return(nerr);
      }
      fprintf(stderr, "PASS    N=%u case word_t!\n\n", N);
      fprintf(stderr, "START free word_t memory\n");
      while (wb)
      {
         wp = wb->next;
         free(wb);
         wb = wp;
      }
      fprintf(stderr, "DONE  free word_t memory\n\n");
      #ifndef WORD_ONLY
         fprintf(stderr, "TESTING N=%u case, i0 item_t\n", N);
         ip = ib;
         ib = mrgSort(ib, gN, g0);
         if (N == 0)
            assert(ib == NULL);
         else if (N == 1)
            assert(ip == ib && ib->next == NULL);
         else
            nerr = tstItems(ib, 0);
         if (nerr)
         {
            fprintf(stderr, "Stopping testing after %u errors with N=%u\n", 
                    nerr, N);
            return(nerr);
         }
         fprintf(stderr, "PASS    N=%u case i0 item_t!\n\n", N);
         fprintf(stderr, "TESTING N=%u case, i1 item_t\n", N);
         ip = ib;
         ib = mrgSort(ib, gN, g1);
         if (N == 0)
            assert(ib == NULL);
         else if (N == 1)
            assert(ip == ib && ib->next == NULL);
         else
            nerr = tstItems(ib, 1);
         if (nerr)
         {
            fprintf(stderr, "Stopping testing after %u errors with N=%u\n", 
                    nerr, N);
            return(nerr);
         }
         fprintf(stderr, "PASS    N=%u case i1 item_t!\n\n", N);
         fprintf(stderr, "START free item_t memory\n");
         while (ib)
         {
            ip = ib->next;
            free(ib);
            ib = ip;
         }
         fprintf(stderr, "DONE  free item_t memory\n\n");
      #endif
   }
   wb = mrgSort(NULL, -1, -1);
   if (wb != NULL)
   {
      fprintf(stderr, "0-length list does not return NULL!\n");
      return(1);
   }
   fprintf(stderr, "\nPASS ALL TESTS!\n\n");
   return(0);
}
