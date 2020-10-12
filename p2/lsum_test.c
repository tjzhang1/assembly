#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
unsigned long long vecsumLL(unsigned int N, unsigned int *V)
{
   unsigned int i;
   unsigned long long sum=0;
   for (i=0; i < N; i++)
      sum += V[i];
   return(sum);
}

int main(int nargs, char **args)
{
   unsigned int N=8000, i, big, last, *X, nerr=0, N32;
   unsigned long long ans, tst;
   unsigned long long lvecsum(unsigned int N, unsigned int *V);
   if (nargs > 1)
      N = atoi(args[1]);
   fprintf(stderr, "TESTING . . . .\n");
   if (lvecsum(0, NULL) != 0)
   {
      fprintf(stderr, "   Can't return 0 for N=0 case!\n");
      nerr++;
   }
/*
 * First, construct a case that will not overflow 32-bits for sanity check
 */
   N32 = (N > 32) ? N : 32;
   X = malloc(sizeof(int)*N32);
   assert(X);
   big = (~0) / N32;
   srand48(N);
   for (i=0; i < N32; i++)
      X[i] = lrand48()%big;
   ans = vecsumLL(N, X);
   tst = lvecsum(N, X);
   if (ans != tst)
   {
     fprintf(stderr, "   32-bit FAIL, GOOD=%llu, BAD=%llu\n", ans, tst);
     nerr++;
   }
   else
      fprintf(stderr, "   32-bit PASS, N=%u, ans=%u\n", N, (unsigned int)ans);
/*
 * Try all cases between 1-24 for cleanup handling
 */
   for (i=1; i <= 24; i++)
   {
      ans = vecsumLL(i, X);
      tst = lvecsum(i, X);
      if (ans != tst)
      {
         fprintf(stderr,"   32-bit N=%u FAIL, GOOD=%llu, BAD=%llu\n",i,ans,tst);
         nerr++;
      }
      else
         fprintf(stderr, "   32-bit N=%u PASS, N=%u, ans=%llu\n", N, i, ans);
   }
/*
 * Next, construct a case that will fail 32-bit accum, but work for 8 accum
 */
   last = (~0) % N;
   if (N >= 8)
      big <<= 3;
   for (i=0; i < last; i++)
      X[i] = big + 1;
   while (i < N)
      X[i++] = big;
   ans = vecsumLL(N, X);
   tst = lvecsum(N, X);
   if (ans != tst)
   {
      fprintf(stderr, "   8x32-bit FAIL, GOOD=%llu, BAD=%llu\n", ans, tst);
      nerr++;
   }
   else
      fprintf(stderr, "   8x32-bit PASS, N=%u, ans=%llu\n", N, ans);
/*
 * Finally, create a case that will fail for 8x32, but pass for 64 bit
 */
   i = big << 3;
   if (i > big) /* may not work if N not large enough */
   {
      big = i;
      for (i=0; i < last; i++)
         X[i] = big + 1;
      while (i < N)
         X[i++] = big;
      ans = vecsumLL(N, X);
      tst = lvecsum(N, X);
      if (ans != tst)
      {
         fprintf(stderr, "   64-bit FAIL, GOOD=%llu, BAD=%llu\n", ans, tst);
         nerr++;
      }
      else
         fprintf(stderr, "   64-bit PASS, N=%u, ans=%llu\n", N, ans);
   }
   free(X);
   if (!nerr)
      printf("ALL TESTS PASS!\n");
   else
      fprintf(stderr, "\nFAILED %u tests!\n\n", nerr);
   return(nerr);
}
