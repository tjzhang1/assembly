/*
float ATL_UDOT(int N, float *X, const int incX, float *Y, const int incY)
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
   return sum0;
}
*/
float ATL_UDOT(int N, float *X, const int incX, float *Y, const int incY)
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
   return sum0;
}
