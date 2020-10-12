#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#if _POSIX_TIMERS > 0
   #include <time.h>
   #define USE_HR 1
#else
   #include <sys/time.h>
   #include <sys/resource.h>
   #define USE_HR 0
#endif
#include "word.h"

double cputime(void) // returns elapsed cputime in seconds as a double
{                    // first call always returns 0, avoid mantissa overflow
   #if USE_HR > 0 // does sys provide high precision timer? */
      struct timespec ts;
   #else
      struct rusage ruse;
   #endif
   static double t0;
   double res;
   static int INIT = 0;

   if (INIT) // common case, where its not first call
   {
      #if USE_HR > 0  
         clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
         res = ts.tv_sec + 1.0e-9*ts.tv_nsec; // translate to double & seconds
      #else
         getrusage(RUSAGE_SELF, &ruse);
         res = ruse.ru_utime.tv_sec + 1.0e-6*ruse.ru_utime.tv_usec;
      #endif
      return(res - t0);                       // subtract frm last call
   }
   #if USE_HR > 0
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&ts);
      t0 = ts.tv_sec + 1.0e-9 * ts.tv_nsec;
   #else
      getrusage(RUSAGE_SELF, &ruse);
      t0 = ruse.ru_utime.tv_sec + 1.0e-6 * ruse.ru_utime.tv_usec;
   #endif
   INIT = 1;
   return(0.0);
}

/*
 * Randomly pivots N-length array, with every element pivoted at least once
 */
void rndPiv(unsigned int N, int *arr)
{
   int i;
   for (i=0; i < N; i++)
   {
      unsigned int k = rand() % N;
      int tmp;
      tmp = arr[i];
      arr[i] = arr[k];
      arr[k] = tmp;
   }
}

int main(int nargs, char **args)
{
   unsigned int i, n;
   word_t *mySort(word_t *bp);

   for (n=1; n < nargs; n++)
   {
      word_t *wp, *wb=NULL;
      unsigned int N;
      double t0, tI, tS, tU;
      char *cp;
      int *idx;

      N = atoi(args[n]);
      /* printf("   INIT N=%u case\n", N); */
      srand(N);
      idx = malloc(sizeof(int)*N);  // get index array
      for (i=0; i < N; i++)         // init them between 0 & N-1
         idx[i] = i;
      rndPiv(N, idx);
      if (N < 30)
      {
         printf("   I=[%u", *idx);
         for (i=1; i < N; i++)
            printf(",%u", idx[i]);
         printf("]\n");
      }
      fflush(stdout);
/*
 *    Allocate N-length lists of words and items, with random int vals
 *    make repeatable by length by seeding with length.
 */
      fflush(stdout);
      t0 = cputime();
/*
 *    Make sure nodes start out randomly placed in allocated space,
 *    otherwise 2nd sort always longer than first due to loss of spacial
 *    locality (perfect before sort!).
 *    NOTE: this illegal, and will only work reliably when 
 *          sizeof(int) == sizeof(ptr)
 */
      cp = calloc(N,sizeof(word_t)); // zero all lens, won't init to 0
      
      for (i=0; i < N; i++)
      {
         int h = idx[i];
         wp = (word_t*)(cp + h*sizeof(word_t));  // not safe
         wp->len = rand();
         wp->next = wb;
         wb = wp;
      }
      free(idx);
      tI = cputime();
      tI = (t0 <= tI) ? tI-t0 : -1.0;
      printf("   INIT   N=%u took %e seconds\n", N, tI);
      printf("   TIMING UNSORTED SORT N=%u case, word_t\n", N);
      fflush(stdout);

      t0 = cputime();
      wb = mySort(wb);
      tU = cputime();
      /* printf("T=%e,%e\n", t0, tU); */
      tU = (t0 <= tU) ? tU-t0 : -1.0;
      printf("   SORT UNSORTED N=%u took %esec, sort/init=%.1f\n", 
             N, tU, tU/tI);
      fflush(stdout);

      t0 = cputime();
      wb = mySort(wb);
      tS = cputime();
      /* printf("T=%e,%e\n", t0, tS); */
      tS = (t0 <= tS) ? tS-t0 : -1.0;
      printf("   SORT   SORTED N=%u took %esec, sort/init=%.1f, "
             "unsort/sort=%.1f\n\n", N, tS, tS/tI, tU/tS);
      free(cp);
   }
   return(0);
}
