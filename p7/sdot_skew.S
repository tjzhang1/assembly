
/*
s0                   r0            r1               r2             r3
float ATL_UDOT(const int N, const float *X, const int incX, const float *Y,
                     SP[0]
                const int incY)
{
   float sum0=0.0,sum1=0.0,sum2=0.0;

   if(N < 6)
   {
      for(;N>0;N--)
      {
         sum0+=X[0]*Y[0];
         X+=1;
         Y+=1;
      }
      return sum0;
   }
   
   float x0=X[0],y0=Y[0],x1=X[1],y1=Y[1],x2=X[2],y2=Y[2];
   float x3=X[3],y3=Y[3],x4=X[4],y4=Y[4],x5=X[5],y5=Y[5];
   X+=6;
   Y+=6;
   for(;N>=12;N-=6)
   {
      sum0+=x0*y0; x0=X[0]; y0=Y[0];
      sum1+=x1*y1; x1=X[1]; y1=Y[1];
      sum2+=x2*y2; x2=X[2]; y2=Y[2];
      sum0+=x3*y3; x3=X[3]; y3=Y[3];
      sum1+=x4*y4; x4=X[4]; y4=Y[4];
      sum2+=x5*y5; x5=X[5]; y5=Y[5];
      X+=6;
      Y+=6;
   }
   sum0+=x0*y0;
   sum1+=x1*y1;
   sum2+=x2*y2;
   sum0+=x3*y3;
   sum1+=x4*y4;
   sum2+=x5*y5;
   X+=6;
   Y+=6;

   sum0+=sum1;
   sum0+=sum2;
   for(N-=6;N>0;N--)
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
#define x0 s3          //value of X[0]
#define y0 s4          //value of Y[0]
#define x1 s5          //value of X[1]
#define y1 s6          //value of Y[1]
#define x2 s7          //value of X[2]
#define y2 s8          //value of Y[2]
#define x3 s9          //value of X[3]
#define y3 s10         //value of Y[3]
#define x4 s11         //value of X[4]
#define y4 s12         //value of Y[4]
#define x5 s13         //value of X[5]
#define y5 s14         //value of Y[5]

.text
.arm
.local LOOP,DONE,LOOP2,CLEANUP1,CLEANUP2
.global ATL_UDOT

ATL_UDOT:
   //current implementation ignores incX and incY
   mov zero,#0         //initialize int reg to 0
   fmsr sum0,zero      //sum0=0.0
   fmsr sum1,zero      //sum1=0.0
   fmsr sum2,zero      //sum2=0.0
   cmp N,#6         //checks there are at least 6 values to accumulate
   bLT CLEANUP1        //if(N<6) skip accumulator loop
   
   //loads 6 values from pX and pY
   flds x0,[pX]        //x0 = *pX;
   flds y0,[pY]        //y0 = *pY;
   flds x1,[pX,#4]     //x1 = pX[1];
   flds y1,[pY,#4]     //y1 = pY[1];
   flds x2,[pX,#8]     //x2 = pX[2];
   flds y2,[pY,#8]     //y2 = pY[2];
   flds x3,[pX,#12]    //x3 = pX[3];
   flds y3,[pY,#12]    //x3 = pY[3];
   flds x4,[pX,#16]    //x4 = pX[4];
   flds y4,[pY,#16]    //y4 = pY[4];
   flds x5,[pX,#20]    //x5 = pX[5];
   flds y5,[pY,#20]    //x5 = pY[5];

   add pX,pX,#24       //pX+=6;
   add pY,pY,#24       //pY+=6;
   cmp N,#12         //checks there are at least 6 more values to load
   bLT CLEANUP2        //if(N<12) skip accumulator loop and sum loaded values
   LOOP:
      fmacs sum0,x0,y0 //sum0 += x0 * y0;
      flds x0,[pX]     //x0 = *pX;
      flds y0,[pY]     //y0 = *pY;

      fmacs sum1,x1,y1 //sum1 += x1 * y1;
      flds x1,[pX,#4]  //x1 = pX[1];
      flds y1,[pY,#4]  //y1 = pY[1];

      fmacs sum2,x2,y2 //sum2 += x2 * y2;
      flds x2,[pX,#8]  //x2 = pX[2];
      flds y2,[pY,#8]  //y2 = pY[2];

      fmacs sum0,x3,y3 //sum0 += x3 * y3;
      flds x3,[pX,#12] //x3 = pX[3];
      flds y3,[pY,#12] //x3 = pY[3];

      fmacs sum1,x4,y4 //sum1 += x4 * y4;
      flds x4,[pX,#16] //x4 = pX[4];
      flds y4,[pY,#16] //y4 = pY[4];

      fmacs sum2,x5,y5 //sum2 += x5 * y5;
      flds x5,[pX,#20] //x5 = pX[5];
      flds y5,[pY,#20] //x5 = pY[5];

      add pX,pX,#24    //pX+=6;
      add pY,pY,#24    //pY+=6;
      subs N,N,#6      //N-=6;
      cmp N,#12
   bGE LOOP            //if(N>=0) continue loop
   
   CLEANUP2:           //still have six loaded values to sum up
   fmacs sum0,x0,y0    //sum0+=x0*y0;
   fmacs sum1,x1,y1    //sum1+=x1*y1;
   fmacs sum2,x2,y2    //sum2+=x2*y2;
   fmacs sum0,x3,y3    //sum0+=x3*y3;
   fmacs sum1,x4,y4    //sum1+=x4*y4;
   fmacs sum2,x5,y5    //sum2+=x5*y5;
   //add pX,pX,#24       //pX+=6;
   //add pY,pY,#24       //pY+=6;

   fadds sum0,sum0,sum1//sum0+=sum1;
   fadds sum0,sum0,sum2//sum0+=sum2;
   subs N,N,#6         //N-=6;
   bLE DONE            //if(N<=0)   skip cleanup

   CLEANUP1:
   LOOP2:              //cleanup loop
      fldmIAs pX!,{x0} //x0 = *pX; also increments pX
      fldmIAs pY!,{y0} //y0 = *pY; also increments pY
      fmacs sum0,x0,y0 //sum0 += x0 * y0;
      subs N,N,#1      //N--;
   bNE LOOP2           //if(N>0) continue loop
   DONE:
bx LR
