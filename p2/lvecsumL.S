#define sum r0                //r0 holds sum
#define pV r1                 //r1 holds pointer to array
#define val r2                //r2 holds values in array
#define N r3                  //r3 holds len and index
#define upper r12             //r12 holds most significant bits
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
   mov sum,#0                 //init lower bits=0
   mov upper,#0               //init upper bits=0
   //check if V is NULL
   cmp pV, #0                 //if(!pV)
   bEQ DONE                   //   don't loop

   //Body
   /*
    //loop to add values to sum, this time with carry
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
      adds sum,sum,val        //add newly lded val to sum, set c if needed
      adc upper,upper,#0      //adds the carry bit to the upper bits 
      adds N,N,#1             //if N isn't 0 (there's still vals remaining)
   bNE LOOP                   //   go back to the start of loop

   DONE: //Epilogue
   mov r1, upper              //move MSB into the return reg
bx lr
