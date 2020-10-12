#include <stdio.h>
#include <stdlib.h>

/*function prototype*/
int irsum(int N, int *pX)
{
   if(N>0)
   {
      if(N>1)
      {
         int midpt;       //holds center of array
         int lsum, rsum;  //vars for left sum and right sum

         midpt = N/2;     //midpt=len of left=starting indx of right
         lsum = irsum(midpt, pX);
         rsum = irsum(N-midpt, pX + midpt);
         return lsum+rsum;
      }
      return pX[0];
   }
   return 0;
}

/*whaley's function*/
/*
int isum(int N, int *X)
{
   if (N > 0)
   {
      int i, isum = *X;
      for (i=1; i < N; i++)
         isum += X[i];
      return(isum);
   }
   return(0);
}

int main(void)
{
   int N=10, *X, i;
   X = malloc(sizeof(int) * N);
   printf("X: [");

   for(i=0; i<N; i++)
   {
      X[i] = i;
      printf("%d, ", i);
   }
   printf("]\n");

   printf("Sum recursive: %d\n", irsum(N, X));
   printf("Sum whaley: %d\n", isum(N, X));
   return 0;
}
*/
