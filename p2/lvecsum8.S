#define sum r0                //r0 holds sum
#define v0 r0                 //r0 also used to hold value 0
#define upper r1              //r1 holds upper bits
#define v1 r1                 //r1 also used to hold value 1
#define v2 r2                 //r2 holds value 2
#define N r3                  //r3 holds len or remaining num of terms
#define S0 r4                 //r4 holds s0
#define S1 r5                 //r5 holds s1
#define S2 r6                 //r6 holds s2
#define S3 r7                 //r7 holds s3
#define S4 r8                 //r8 holds s4
#define S5 r9                 //r9 holds s5
#define S6 r10                //r10 holds s6
#define S7 r11                //r11 holds s7
#define pV r12                //r12 holds ptr to V
#define v3 r14                //r14 holds value 3
/*
r0,r1                                  r0               r1
unsigned long long lvecsum(unsigned int N, unsigned int *V)
{
   unsigned long long sum=0;
   unsigned int s0=0, s1=0, s2=0, s3=0, s4=0, s5=0, s6=0, s7=0;
   for (N=8-N; N<=0; N+=8)
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
   for (N-=8; N<0; N++)
   {
      sum += *(V++);
   }
   sum += s0+s1+s2+s3+s4+s5+s6+s7;
   return(sum);
}
*/

.text
.arm
.local LOOP1,LOOP2,DONE
.global lvecsum

lvecsum:
   //Prologue
   mov N,r0                   //move N so sum can be stored in r0
   mov pV,r1                  //move pV so upper bits can be stored in r1
   mov sum,#0                 //init lower bits=0
   mov upper,#0               //init upper bits=0
   //check if V is NULL
   cmp pV, #0                 //if(!pV)
   bEQ DONE                   //   don't loop
   //save caller's r4-r11 and r14, and save values for sum and upper
   stmDB SP!,{sum,upper,r4-r11,r14} //makes space for s0-s7 and v0-v3

   //Body
   //initialize partial sums to 0
   mov S0,#0                  //S0=0
   mov S1,#0                  //S1=0
   mov S2,#0                  //S2=0
   mov S3,#0                  //S3=0
   mov S4,#0                  //S4=0
   mov S5,#0                  //S5=0
   mov S6,#0                  //S6=0
   mov S7,#0                  //S7=0
   /*
    //this loop updates the partial sums 8 at a time
    //N starts negative and acts as counter for remaining num of terms
    //start at 8-N to avoid out-of-bounds loads if N is not multiple of 8
      for (N=8-N; N<=0; N+=8)
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
   */
   rsbs N,N,#8                //N=8-N; if <=0 then enter loop
   LOOP1:
      ldmLEIA pV!,{v0,v1,v2,v3}  //load next four values from the pointer
      addLE S0,S0,v0          //add next val to s0
      addLE S1,S1,v1          //add next val to s1
      addLE S2,S2,v2          //add next val to s2
      addLE S3,S3,v3          //add next val to s3
      ldmLEIA pV!,{v0,v1,v2,v3}  //load next four values from the pointer
      addLE S4,S4,v0          //add next val to s4
      addLE S5,S5,v1          //add next val to s5
      addLE S6,S6,v2          //add next val to s6
      addLE S7,S7,v3          //add next val to s7
      addLEs N,N,#8           //N+=8; if N<=0,
   bLE LOOP1                  //         then continue loop
   //now v0-v3 (which are r0-r2,r14) are free for use
   ldmIA SP!,{sum,upper}      //restore sum and upper regs
   /*
    //loop to add remaining vals to sum; there will always be fewer than 8 left
    //N is in range [1,8]; sub 8 to find number of remaining vals
      for (N-=8; N<0; N++)
      {
         sum += *(V++);
      }
   */
   subs N,N,#8                //sub 8 from N to get num of remaining vals
   //if N is already zero, don't loop
   LOOP2:
      ldmNEIA pV!,{v2}        //load next value of pV into v2
      addNE sum,sum,v2        //sum+=v2
      addNEs N,N,#1           //N++
   bNE LOOP2                  //branch until hits zero
   //add up partial sums
   adds sum,sum,S0            //sum+=s0, sets c if needed
   adc upper,upper,#0         //adds the carry bit to the upper bits
   adds sum,sum,S1            //sum+=s1, sets c if needed
   adc upper,upper,#0         //adds the carry bit to the upper bits
   adds sum,sum,S2            //sum+=s2, sets c if needed
   adc upper,upper,#0         //adds the carry bit to the upper bits
   adds sum,sum,S3            //sum+=s3, sets c if needed
   adc upper,upper,#0         //adds the carry bit to the upper bits
   adds sum,sum,S4            //sum+=s4, sets c if needed
   adc upper,upper,#0         //adds the carry bit to the upper bits
   adds sum,sum,S5            //sum+=s5, sets c if needed
   adc upper,upper,#0         //adds the carry bit to the upper bits
   adds sum,sum,S6            //sum+=s6, sets c if needed
   adc upper,upper,#0         //adds the carry bit to the upper bits
   adds sum,sum,S7            //sum+=s7, sets c if needed
   adc upper,upper,#0         //adds the carry bit to the upper bits

   //Epilogue
   ldmIA SP!,{r4-r11,r14}     //restore caller's regs and dealloc stack space
   DONE:
bx lr
