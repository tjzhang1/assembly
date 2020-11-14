void dgemm_asg(int ta, int tb, int M, int N, int K, double alpha, double *A,
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
