#include <stdio.h>
#include <stdlib.h>

void printPtr(int *ptr, int len);
int *init(int len);

int main(void)
{
   int N=10;
   int s0=0,s1=0,s2=0,s3=0,s4=0,s5=0,s6=0,s7=0;
   int sum=0;
   int *V = init(N), *p;
   p = V;

   printPtr(V, N);
   for (N=N-8; N>=0; N-=8)
   {
      s0 += *(V++);
      s1 += *(V++);
      s2 += *(V++);
      s3 += *(V++);
      s4 += *(V++);
      s5 += *(V++);
      s6 += *(V++);
      s7 += *(V++);
   }
   for (N+=8; N>0; N--)
   {
      sum += *(V++);
   }
   sum += s0+s1+s2+s3+s4+s5+s6+s7;
   printf("Sum=%d\n", sum);
   free(p);
   return 0;
}

void printPtr(int *ptr, int len)
{
   int i;
   for(i=0;i<len;i++)
   {
      printf("%d, ", ptr[i]);
   }
   printf("\n");
}

int *init(int len)
{
   int i;
   int *p;

   p = malloc(sizeof(int)*len);
   for(i=0;i<len;i++)
   {
      p[i]=i+1;
   }
   return p;
}
