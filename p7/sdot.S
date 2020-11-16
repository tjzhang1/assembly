
/*
s0                   r0            r1               r2             r3
float ATL_UDOT(const int N, const float *X, const int incX, const float *Y,
                     SP[0]
                const int incY)
{
   float sum0=0.0,sum1=0.0,sum2=0.0;

   for(N-=6;N>=0;N-=6)
   {
      sum0+=X[0]*Y[0];
      sum1+=X[1]*Y[1];
      sum2+=X[2]*Y[2];
      sum0+=X[3]*Y[3];
      sum1+=X[4]*Y[4];
      sum2+=X[5]*Y[5];
      X+=6;
      Y+=6;
   }
   sum0+=sum1;
   sum0+=sum2;
   for(N+=6;N>0;N--)
   {
      sum0 += X[0]*Y[0];
      X+=1;
      Y+=1;
   }
   return sum;
}
*/

//ASSUMES: *X and *Y are arrays length N, N>=1

#define N r0           //holds num of remaining elements
#define pX r1          //ptr to first element in X array
#define zero r2        //this is zero register, to init accumulators
#define pY r3          //ptr to first element in Y array
#define sum0 s0        //accumulator 0 
#define sum1 s1        //accumulator 1
#define sum2 s2        //accumulator 2
#define xV s3          //value of element from X array
#define yV s4          //value of element from Y array

.text
.arm
.local LOOP,DONE,LOOP2,CLEANUP
.global ATL_UDOT

ATL_UDOT:
   //current implementation ignores incX and incY
   mov zero,#0         //initialize int reg to 0
   fmsr sum0,zero      //sum0=0.0
   fmsr sum1,zero      //sum1=0.0
   fmsr sum2,zero      //sum2=0.0
   subs N,N,#6         //checks there are at least 6 values to accumulate
   bLT CLEANUP         //if(N<0) skip accumulator loop

   /*
   *for(N-=6;N>=0;N-=6)
   *{
   *   sum0+=X[0]*Y[0];
   *   sum1+=X[1]*Y[1];
   *   sum2+=X[2]*Y[2];
   *   sum0+=X[3]*Y[3];
   *   sum1+=X[4]*Y[4];
   *   sum2+=X[5]*Y[5];
   *   X+=6;
   *   Y+=6;
   *}
   */ 
   LOOP:
      fldmIAs pX!,{xV} //xV = *pX; also increments pX
      fldmIAs pY!,{yV} //yV = *pY; also increments pY
      fmacs sum0,xV,yV //sum0 += xV * yV;

      fldmIAs pX!,{xV} //xV = *pX; also increments pX
      fldmIAs pY!,{yV} //yV = *pY; also increments pY
      fmacs sum1,xV,yV //sum1 += xV * yV;

      fldmIAs pX!,{xV} //xV = *pX; also increments pX
      fldmIAs pY!,{yV} //yV = *pY; also increments pY
      fmacs sum2,xV,yV //sum2 += xV * yV;

      fldmIAs pX!,{xV} //xV = *pX; also increments pX
      fldmIAs pY!,{yV} //yV = *pY; also increments pY
      fmacs sum0,xV,yV //sum0 += xV * yV;

      fldmIAs pX!,{xV} //xV = *pX; also increments pX
      fldmIAs pY!,{yV} //yV = *pY; also increments pY
      fmacs sum1,xV,yV //sum1 += xV * yV;

      fldmIAs pX!,{xV} //xV = *pX; also increments pX
      fldmIAs pY!,{yV} //yV = *pY; also increments pY
      fmacs sum2,xV,yV //sum2 += xV * yV;

      subs N,N,#6      //N-=6;
   bGT LOOP            //if(N>=0) continue loop

   fadds sum0,sum0,sum1//sum0+=sum1;
   fadds sum0,sum0,sum2//sum0+=sum2;

   CLEANUP:
   adds N,N,#6         //N+=6; if(N<=0)
   bLE DONE            //   skip cleanup
   /*
   *for(N+=6;N>0;N--)
   *{
   *   sum0 += X[0]*Y[0];
   *   X+=1;
   *   Y+=1;
   *}
   */
   LOOP2:              //cleanup loop
      fldmIAs pX!,{xV} //xV = *pX; also increments pX
      fldmIAs pY!,{yV} //yV = *pY; also increments pY
      fmacs sum0,xV,yV //sum0 += xV * yV;
      subs N,N,#1      //N--;
   bNE LOOP2           //if(N>0) continue loop
   DONE:
bx LR
