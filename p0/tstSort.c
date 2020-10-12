#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "word.h"

unsigned int cntNodes(word_t *wb)
{
   unsigned int cnt=0;
   while (wb)
   {
      cnt++;
      wb = wb->next;
   }
   return(cnt);
}

int main(int nargs, char **args)
{
   unsigned int i, n;
   word_t *wp, *wb=NULL;
   word_t *mySort(word_t *lb);


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
         #if 0
            wp->len = ip->i0 = ip->i1 = N-i;
         #else
            wp->len = rand();
         #endif
      }
      wp = wb;
      wb = mySort(wb);
      assert(cntNodes(wb) == N);
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
   }
   wb = mySort(NULL);
   if (wb != NULL)
   {
      fprintf(stderr, "0-length list does not return NULL!\n");
      return(1);
   }
   fprintf(stderr, "\nPASS ALL TESTS!\n\n");
   return(0);
}
