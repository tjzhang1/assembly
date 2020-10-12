#define sum r0                //r0 holds sum
#define pV r1                 //r1 holds pointer to array
#define val r2                //r2 holds values in array
#define N r3                  //r3 holds len and index
/*
r0,r1                                  r0               r1
unsigned long long lvecsum(unsigned int N, unsigned int *V)
{
   unsigned int i;
   unsigned long long sum=0;
   for (i=0; i < N; i++)
      sum += V[i];
   return(sum);
}
*/

.text
.arm
.local LOOP, DONE
.global lvecsum

lvecsum:
   //Prologue
   mov N,r0                   //move N so sum can be stored in r0
   mov sum,#0                 //init sum=0
   //check if V is NULL
   cmp pV, #0                 //if(!pV)
   bEQ DONE                   //   don't loop

   //Body
   /*
    //loop to add values to sum
      V+=N;
      for (N=-N; N<0; N++)
         sum += V[N];
   */
   //move pointer to end of array
   add pV,pV,N,LSL #2         //V+=N
   //set index to -N
   rsb N,N,#0                 //N=0-N=-N
   LOOP:
      ldr val,[pV,N,LSL #2]   //val=V[N]
      add sum,sum,val         //add newly loaded value to sum
      adds N,N,#1             //if N isn't zero (there's still nums remaining)
   bNE LOOP                   //   then go back to start of loop

   DONE: //Epilogue
   mov r1,#0                  //clear high bits
bx lr
