/*
                   r0      r1     r2     r3  SP[0]         d0          SP[4]
void dgemm_asg(int ta, int tb, int M, int N, int K, double alpha, double *A,
                 SP[8]     SP[12]   SP[16]         d1       SP[20]    SP[24]
               int lda, double *B, int ldb, double beta, double *C, int ldc)
{
   //assume A is transposed = AT
   int i,j,k;

   for (i=0; i<N; i++) //loop through cols of B
   {
      for (j=0; j<M; j++) //loop through cols of AT (rows of A)
      {
         double valC, product=0.0;

         if(beta == 0.0)
            valC = beta;
         else
            valC = beta*C[i*ldc+j];
         //dot product
         for(k=0; k<K; k+=1)  //loop through elements in A/B arrays for dot
         {
            product += A[j*lda+k]*B[i*ldb+k];
         }
         C[i*ldc+j] = (alpha*product+valC); //C col major
      }
   }
}
*/

.arm
.text
.local LOOPB,LOOPA,LOOPK
.global dgemm_asg

#define K r0              //shared length between A,B
#define lda r1            //stride length of A, keeps track of current col A
#define M r2              //num cols of At (= num rows of A)
#define N r3              //num cols of B
#define ldb r4            //stride length of B, keeps track of current col B
#define ldc r5            //stride length of C, keeps track of current col C
#define pA r6             //ptr to A
#define pB r7             //ptr to B
#define pC r8             //ptr to C
#define i r9              //index through N (cols of B)
#define j r10             //index through M (cols of At)
#define k r12             //index through K (length of arrays)
#define alpha d0
#define beta d1
#define vA d2             //value of A
#define vB d3             //value of B
#define vC d4             //value of C
#define dp d5             //dot product

#define FSZ 7*4           //saving four int regs
#define argK FSZ          //K = SP[0]
#define argPA FSZ+4       //pA = SP[4]
#define argLda FSZ+8      //lda = SP[8]
#define argPB FSZ+12      //pB = SP[12]
#define argLdb FSZ+16     //ldb = SP[16]
#define argPC FSZ+20      //pC = SP[20]
#define argLdc FSZ+24     //ldc = SP[24]

dgemm_asg:
   stmDB SP!,{r4-r10}      //save caller's regs
   ldr K,[SP,#argK]        //load K
   ldr pA,[SP,#argPA]      //load pA
   ldr lda,[SP,#argLda]    //load lda
   ldr pB,[SP,#argPB]      //load pB
   ldr ldb,[SP,#argLdb]    //load ldb
   ldr pC,[SP,#argPC]      //load pC
   ldr ldc,[SP,#argLdc]    //load ldc
   
   //Assume N,M,K >= 1
   mov i,#0                //init i=0;
   LOOPB:                  //loop through N cols of B
   //{
      mov j,#0             //init j=0;
      LOOPA:                  //loop through M cols of At
      //{
         mov k,#0             //init k=0;
         fmdlr dp,k           //mov 0 to low bits of dp
         fmdhr dp,k           //init product = 0.0
         fcmpd beta,dp        //vC = (beta == 0.0) ? 0.0:beta*C[0]
         fmstat               //   set iCC
         fcpydEQ vC,beta      //   vC = 0.0 if beta == 0.0
         flddNE vC,[pC]       //   vC = *C if beta != 0.0
         fmuldNE vC,vC,beta   //   vC = beta*C[0] if beta != 0.0
         LOOPK:                  //loop through K elements for dot
         //{
	    //using ptr arithmetic for indexing b/c can't use fldr with iregs
            fldmIAd pA!,{vA}     //load A and increment ptr
            fldmIAd pB!,{vB}     //load B and increment ptr
            fmacd dp,vA,vB       //dp += vA*vB
            add k,k,#1           //k++;
            cmp k,K              //if(k<K)
         bLT LOOPK               //   continue looping through At/B arrays
         //}
         sub pA,pA,k,LSL #3   //A-=k, removes accumulated offset from LOOPK
         sub pB,pB,k,LSL #3   //B-=k, removes accumulated offset from LOOPK
         fmacd vC,dp,alpha    //vC += alpha*dp;
         fstd vC,[pC]         //C[i*ldc+j] = alpha*dp+beta;
         add pA,pA,lda,LSL #3 //increments A by one stride
         add pC,pC,#8         //increments C
         add j,j,#1           //j++;
         cmp j,M              //if(j<M)
      bLT LOOPA               //   continue looping through At
      //}
      sub pC,pC,j,LSL #3   //C-=j, removes accumulated offset from LOOPA
      mul j,j,lda          //j *= lda
      sub pA,pA,j,LSL #3   //A-=j*lda, removes accumulated offset from LOOPA
      add pB,pB,ldb,LSL #3 //increments pB by one stride
      add pC,pC,ldc,LSL #3 //increments pC by one stride
      add i,i,#1           //i++;
      cmp i,N              //if(i<N)
   bLT LOOPB               //   continue looping through B
   //}
   ldmIA SP!,{r4-r10}      //restore caller's regs
bx LR
