
/*
d0                   r0            r1               r2               r3
double ATL_UDOT(const int N, const double *X, const int incX, const double *Y,
                     SP[0]
                const int incY)
{
   double sum = 0.0;
   int i;
   
   for(i=0;i<N;i++)
      sum+=X[i]*Y[i];
   return sum;
}
*/

//ASSUMES: *X and *Y are arrays length N, N>=1

#define N r0           //holds num of remaining elements
#define pX r1          //ptr to first element in X array
#define pY r3          //ptr to first element in Y array
#define sum d0         //accumulator of type double
#define xV d1          //value of element from X array
#define yV d2          //value of element from Y array

.text
.arm
.local LOOP,DONE
.global ATL_UDOT

ATL_UDOT:
   //current implementation ignores incX and incY
   fldmIAd pX!,{xV}    //xV = *pX; also increments pX
   fldmIAd pY!,{yV}    //yV = *pY; also increments pY
   fmuld sum,xV,yV     //sum = xV * yV;
   subs N,N,#1         //N--;
   bEQ DONE            //if(!N) skip loop

   LOOP:
      fldmIAd pX!,{xV} //xV = *pX; also increments pX
      fldmIAd pY!,{yV} //yV = *pY; also increments pY
      fmacd sum,xV,yV  //sum += xV * yV;
      subs N,N,#1      //N--;
   bNE LOOP            //if(N) continue loop
   DONE:
bx LR
