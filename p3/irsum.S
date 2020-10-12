
/*
r0           r0       r1
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
*/

.text
.arm
.global irsum

#define mdpt r4
#define lsum r5
#define N r6
#define pX r7

irsum:
   //Prologue
   cmp r0,#1                  //if(N<=1) return some value specified below:
   movLT r0,#0                //if(N<1) return 0
   ldrEQ r0,[r1]              //if(N==1) return pX[0]
   bxLE lr                    //return value specified above
   
   stmDB SP!,{r4-r7,LR}       //store caller's regs
   mov N,r0                   //move N to callee saved reg
   mov pX,r1                  //move pX to callee saved reg

   //Body
   mov mdpt,N,LSR #1          //mdpt=N/2, save to reg for use after recursion
   //left call
   mov r0,mdpt                //arg0=mdpt, which is len of left half
   //pX is still in r1 at this point
   bl irsum                   //irsum(mdpt,pX)
   mov lsum,r0                //save returned left sum
   //right call
   sub r0,N,mdpt              //arg0=N-mdpt, which is len of right half
   add r1,pX,mdpt,LSL #2      //arg1=pX+mdpt, starting point for right half
   bl irsum                   //irsum(N-mdpt,pX+mdpt)
   add r0,r0,lsum             //add left sum into right sum
   
   //Epilogue
   ldmIA SP!, {r4-r7,LR}      //restore caller's regs
bx lr

