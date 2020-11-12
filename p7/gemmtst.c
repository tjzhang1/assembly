#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define time00 ATL_cputime
#ifndef TYPE
   #define EPS 1.0e-15
   #define TYPE double
   #define PRE d
   #define UPR d
   #define PREU D
   #define PATL ATL_d
   #define PATU ATLU_d
   #define UATL ATLU_d
   #define CBLA cblas_d
   #define PATLU ATL_d
   #define ATL_rone   1.0
   #define ATL_rnone -1.0
   #define ATL_rzero   0.0
   #define ATL_typify(m_) m_
   #define SHIFT
   #define SCALAR TYPE
   #define SADD &
   #define SVAL
   #define SVVAL *
   #define SCALAR_IS_ONE(M_scalar) ((M_scalar) == ATL_rone)
   #define SCALAR_IS_NONE(M_scalar) ((M_scalar) == ATL_rnone)
   #define SCALAR_IS_ZERO(M_scalar) ((M_scalar) == ATL_rzero)
   #define TREAL
   #define ATL_sizeof 8
   #define ATL_MulBySize(i_) ((i_)<<3)
#endif
#define Mjoin(pre, nam) my_join(pre, nam)
#define my_join(pre, nam) pre ## nam
#define Mstr2(m) # m
#define Mstr(m) Mstr2(m)
#define Mabs(x) ( (x) >= 0 ? (x) : -(x) )
#define Mmax(x, y) ( (x) > (y) ? (x) : (y) )
#define Mmin(x, y) ( (x) > (y) ? (y) : (x) )
#define Mlowcase(C) ( ((C) > 64 && (C) < 91) ? (C) | 32 : (C) )
#define Mupcase(C) ( ((C) > 96 && (C) < 123) ? (C) & 0xDF : (C) )
#define ATL_assert(n_) \
{ \
   if (!(n_)) \
   { \
      fprintf(stderr, "assertion %s failed, line %d of file %s\n", \
                 Mstr(n_), __LINE__, __FILE__); \
   } \
}
#define ATL_Cachelen 32
#define ATL_MulByCachelen(N_) ( (N_) << 5 )
#define ATL_DivByCachelen(N_) ( (N_) >> 5 )
#define ATL_AlignPtr(vp) \
   (void*) (ATL_Cachelen + ATL_MulByCachelen(ATL_DivByCachelen((size_t) (vp))))

#define dumb_seed(iseed_) srand(iseed_)
#ifndef RAND_MAX  /* rather dangerous non-ansi workaround */
   #define RAND_MAX ((unsigned long)(1<<30))
#endif
#define dumb_rand() ( 0.5 - ((double)rand())/((double)RAND_MAX) )



enum ATLAS_TRANS {AtlasNoTrans=111, AtlasTrans=112, AtlasConjTrans=113,
                  AtlasConj=114};

#ifndef L2SIZE
   #define L2SIZE 4194304
#endif

#define dgemmRCW dgemm_tst2
void dgemm
   (int ta, int tb, int M, int N, int K, double alpha, const double *A,
    int lda, const double *B, int ldb, double beta, double *C, int ldc);
void dgemm_asg
   (int ta, int tb, int M, int N, int K, double alpha, const double *A,
    int lda, const double *B, int ldb, double beta, double *C, int ldc);

static void dgemm_tst2
   (int ta, int tb, int M, int N, int K, double alpha, const double *A,
    int lda, const double *B, int ldb, double beta, double *C, int ldc)
{
   int i, j;

   if (M < 1 || N < 1)
      return;
   if (alpha == 0.0 || K < 1)
   {
      if (beta != 1.0)
      {
         if (beta == 0.0)
         {
            for (j=0; j < N; j++)
            {
               for (i=0; i < M; i++)
                  C[i+j*ldc] = 0.0;
            }
         }
         else
         {
            for (j=0; j < N; j++)
            {
               for (i=0; i < M; i++)
                  C[i+j*ldc] *= beta;
            }
         }
      }
      return;
   }
   dgemm(ta, tb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

#if defined(UseClock)
   #include <time.h>
   double ATL_cputime(void)
   {
      clock_t t1;
      static int INIT=0;
      static clock_t t0;
      static const double CPS = 1.0 / (1.0*CLOCKS_PER_SEC);
      double d;

      if (INIT)
      {
         t1 = clock() - t0;
         d = t1 * CPS;
         return(d);
      }
      INIT = 1;
      t0 = clock();
      return(0.0);
   }
#elif defined(UseTimes)
   #include <sys/times.h>
   #include <unistd.h>
   double ATL_cputime(void)
   {
      struct tms ts;
      static double ClockTick=0.0;

      if (ClockTick == 0.0) ClockTick = 1.0 / ((double) sysconf(_SC_CLK_TCK));
      times(&ts);
      return( ((double) ts.tms_utime) * ClockTick );
   }
#elif defined(SUN_HR) /* use sun high resolution timers */
   #include <sys/time.h>
   double ATL_cputime(void)
   {
      return(gethrvtime()*1.0e-9);
   }
#else
   #include <sys/time.h>
   #include <sys/resource.h>
   double ATL_cputime(void)
   {
      struct rusage ruse;
      getrusage(RUSAGE_SELF, &ruse);
      return( (double)(ruse.ru_utime.tv_sec+ruse.ru_utime.tv_usec/1000000.0) );
   }
#endif


#ifdef ATL_DeclareSlens
F77_INTEGER ATL_Slen1, ATL_Slen2;
#endif
double time00();


#include <stdarg.h>
void ATL_xerbla(int p, char *rout, char *form, ...)
{
   va_list argptr;

   va_start(argptr, form);
#ifdef GCCWIN
   if (p) printf("Parameter %d to routine %s was incorrect\n", p, rout);
   vprintf(form, argptr);
#else
   if (p)
      fprintf(stderr, "Parameter %d to routine %s was incorrect\n", p, rout);
   vfprintf(stderr, form, argptr);
#endif
   va_end(argptr);
   exit(-1);
}

void Mjoin(PATL,zero)(const int N, TYPE *X, const int incX)
/*
 * X <- 0
 */
{
   int i, n;
   if (incX == 1)
   {
      n = N SHIFT;
      i = n >> 5;
      if (i)
      {
         n -= (i << 5);
         do
         {
            *X = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = X[7] = X[8] = X[9] =
            X[10] = X[11] = X[12] = X[13] = X[14] = X[15] = X[16] = X[17] =
            X[18] = X[19] = X[20] = X[21] = X[22] = X[23] = X[24] = X[25] =
            X[26] = X[27] = X[28] = X[29] = X[30] = X[31] = ATL_rzero;
            X += 32;
         }
         while(--i);
      }
      if (n >> 4) /* >= 16 */
      {
         *X = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = X[7] = X[8] = X[9] =
         X[10] = X[11] = X[12] = X[13] = X[14] = X[15] = ATL_rzero;
         X += 16;
         n -= 16;
      }
      if (n >> 3) /* >= 8 */
      {
         *X = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = X[7] = ATL_rzero;
         X += 8;
         n -= 8;
      }
      switch(n)
      {
         case 1:
            *X = ATL_rzero;
            break;
         case 2:
            *X = X[1] = ATL_rzero;
            break;
         case 3:
            *X = X[1] = X[2] = ATL_rzero;
            break;
         case 4:
            *X = X[1] = X[2] = X[3] = ATL_rzero;
            break;
         case 5:
            *X = X[1] = X[2] = X[3] = X[4] = ATL_rzero;
            break;
         case 6:
            *X = X[1] = X[2] = X[3] = X[4] = X[5] = ATL_rzero;
            break;
         case 7:
            *X = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = ATL_rzero;
            break;
         default:;
      }
   }
   else
   {
      #ifdef TREAL
         for (i=N; i; i--, X += incX) *X = ATL_rzero;
      #else
         for (n=incX<<1, i=N; i; i--, X += n) *X = X[1] = ATL_rzero;
      #endif
   }
}
double ATL_flushcache(int size)
/*
 * flush cache by reading enough mem; note that if the compiler gets
 * really smart, may be necessary to make vp a global variable so it
 * can't figure out it's not being modified other than during setup;
 * the fact that ATL_dzero is external will confuse most compilers
 */
{
  static void *vp=NULL;
  static int N = 0;
  double *cache;
  double dret=0.0;
  int i;

  if (size < 0) /* flush cache */
  {
     ATL_assert(vp);
     cache = ATL_AlignPtr(vp);
     if (N > 0) for (i=0; i != N; i++) dret += cache[i];
  }
  else if (size > 0) /* initialize */
  {
     vp = malloc(ATL_Cachelen + size);
     ATL_assert(vp);
     N = size / sizeof(double);
     cache = ATL_AlignPtr(vp);
     ATL_dzero(N, cache, 1);
  }
  else if (size == 0) /* free cache */
  {
     if (vp) free(vp);
     vp = NULL;
     N = 0;
  }
  return(dret);
}
static void ATL_ladd
(
   int *                      J,
   int *                      K,
   int *                      I
)
{
/*
 * Purpose
 * =======
 *
 * ATL_ladd adds  without carry two long positive integers  K and J  an
 * put the result into I.  The long integers  I, J, K are encoded on 31
 * bits using an array of 2 integers.  The 16-lower bits  are stored  i
 * the  first  entry  of each array,  the 15-higher bits  in the second
 * entry.
 *
 * Arguments
 * =========
 *
 * J       (local input)                 int *
 *         On entry, J is an integer array of dimension 2 containing the
 *         encoded long integer J.
 *
 * K       (local input)                 int *
 *         On entry, K is an integer array of dimension 2 containing the
 *         encoded long integer K.
 *
 * I       (local output)                int *
 *         On entry, I is an integer array of dimension 2. On exit, this
 *         array contains the encoded long integer result.
 *
 * ---------------------------------------------------------------------
 */
   int                        itmp0 = K[0] + J[0], itmp1;
/*
 *    K[1] K[0] K  I[0]  = (K[0]+J[0]) % 2^16
 *    0XXX XXXX    carry = (K[0]+J[0]) / 2^16
 *
 * +  J[1] J[0] J  I[1] = K[1] + J[1] + carry
 *    0XXX XXXX    I[1] = I[1] % 2^15
 *    -------------
 *    I[1] I[0]
 *    0XXX XXXX I
 */
   itmp1 = itmp0 >> 16;         I[0] = itmp0 - ( itmp1 << 16 );
   itmp0 = itmp1 + K[1] + J[1]; I[1] = itmp0 - (( itmp0 >> 15 ) << 15);
}

static void ATL_lmul
(
   int *                      K,
   int *                      J,
   int *                      I
)
{
/*
 * Purpose
 * =======
 *
 * ATL_lmul multiplies  without carry two long positive integers K and J
 * and put the result into I.  The long integers  I, J, K are encoded on
 * 31 bits using an array of 2 integers. The 16-lower bits are stored in
 * the first entry of each array, the 15-higher bits in the second entry
 * of each array. For efficiency purposes, the  intrisic modulo function
 * is inlined.
 *
 * Arguments
 * =========
 *
 * K       (local input)                 int *
 *         On entry, K is an integer array of dimension 2 containing the
 *         encoded long integer K.
 *
 * J       (local input)                 int *
 *         On entry, J is an integer array of dimension 2 containing the
 *         encoded long integer J.
 *
 * I       (local output)                int *
 *         On entry, I is an integer array of dimension 2. On exit, this
 *         array contains the encoded long integer result.
 *
 * ---------------------------------------------------------------------
 */
   static int                 ipow30 = ( 1 << 30 );
   int                        kt, lt;
/*
 *    K[1] K[0] K  kt = K[0]*J[0]
 *    0XXX XXXX    if(kt < 0) kt += 2^31
 * x               I[0] = kt % 2^16
 *                 lt = K[0]*J[1] + K[1]*J[0]
 *    J[1] J[0] J  if(lt < 0) lt += 2^31
 *    0XXX XXXX    kt = (kt / 2^16) + lt
 * --------------  if(kt < 0) kt += 2^31
 *    I[1] I[0]    I[1] = kt % 2^15
 *    0XXX XXXX I
 */
   kt   = K[0] * J[0]; if( kt < 0 ) kt = ( kt + ipow30 ) + ipow30;
   I[0] = kt - ( ( kt >> 16 ) << 16 );
   lt   = K[0] * J[1] + K[1] * J[0];
   if( lt < 0 ) lt = ( lt + ipow30 ) + ipow30;
   kt = ( kt >> 16 ) + lt;
   if( kt < 0 ) kt = ( kt + ipow30 ) + ipow30;
   I[1] = kt - ( ( kt >> 15 ) << 15 );
}

static void ATL_setran
(
   const int                  OPTION,
   int *                      IRAN
)
{
/*
 * Purpose
 * =======
 *
 * ATL_setran initializes  the random generator with the encoding of the
 * first number X(0) in the sequence,  and the constants a and c used to
 * compute the next element in the sequence: X(n+1) = a*X(n) + c.  X(0),
 * a and c are stored in the static variables  irand, ias and ics.  When
 * OPTION is 0 (resp. 1 and 2),  irand  (resp. ia and ic)  is set to the
 * values of the input array IRAN.  When OPTION is 3, IRAN is set to the
 * current value of irand, and irand is then incremented.
 *
 * Arguments
 * =========
 *
 * OPTION  (local input)                 const int
 *         On entry, OPTION  is an integer that specifies the operations
 *         to be performed on the random generator as specified above.
 *
 * IRAN    (local input/output)          int *
 *         On entry,  IRAN is an array of dimension 2, that contains the
 *         16-lower and 15-higher bits of a random number.
 *
 * ---------------------------------------------------------------------
 */
   static int                 ias[2], ics[2], irand[2];
   int                        j[2];

   if(      OPTION == 3 )
   {                                       /* return current value */
      IRAN[0] = irand[0]; IRAN[1] = irand[1];
      ATL_lmul( irand, ias, j );         /* j     = irand * ias;   */
      ATL_ladd( j, ics, irand );         /* irand = j     + ics;   */
   }
   else if( OPTION == 0 ) { irand[0] = IRAN[0]; irand[1] = IRAN[1]; }
   else if( OPTION == 1 ) { ias  [0] = IRAN[0]; ias  [1] = IRAN[1]; }
   else if( OPTION == 2 ) { ics  [0] = IRAN[0]; ics  [1] = IRAN[1]; }
}

static void ATL_xjumpm
(
   const int                  JUMPM,
   int *                      MULT,
   int *                      IADD,
   int *                      IRANN,
   int *                      IRANM,
   int *                      IAM,
   int *                      ICM
)
{
/*
 * Purpose
 * =======
 *
 * ATL_xjumpm computes  the constants  A and C  to jump JUMPM numbers in
 * the random sequence: X(n+JUMPM) = A*X(n)+C.  The constants encoded in
 * MULT and IADD  specify  how to jump from one entry in the sequence to
 * the next.
 *
 * Arguments
 * =========
 *
 * JUMPM   (local input)                 const int
 *         On entry,  JUMPM  specifies  the  number  of entries  in  the
 *         sequence to jump over. When JUMPM is less or equal than zero,
 *         A and C are not computed, IRANM is set to IRANN corresponding
 *         to a jump of size zero.
 *
 * MULT    (local input)                 int *
 *         On entry, MULT is an array of dimension 2,  that contains the
 *         16-lower  and 15-higher bits of the constant  a  to jump from
 *         X(n) to X(n+1) = a*X(n) + c in the random sequence.
 *
 * IADD    (local input)                 int *
 *         On entry, IADD is an array of dimension 2,  that contains the
 *         16-lower  and 15-higher bits of the constant  c  to jump from
 *         X(n) to X(n+1) = a*X(n) + c in the random sequence.
 *
 * IRANN   (local input)                 int *
 *         On entry, IRANN is an array of dimension 2. that contains the
 *         16-lower and 15-higher bits of the encoding of X(n).
 *
 * IRANM   (local output)                int *
 *         On entry,  IRANM  is an array of dimension 2.   On exit, this
 *         array  contains respectively  the 16-lower and 15-higher bits
 *         of the encoding of X(n+JUMPM).
 *
 * IAM     (local output)                int *
 *         On entry, IAM is an array of dimension 2. On exit, when JUMPM
 *         is  greater  than  zero,  this  array  contains  the  encoded
 *         constant  A  to jump from  X(n) to  X(n+JUMPM)  in the random
 *         sequence. IAM(0:1)  contains  respectively  the  16-lower and
 *         15-higher  bits  of this constant  A. When  JUMPM  is less or
 *         equal than zero, this array is not referenced.
 *
 * ICM     (local output)                int *
 *         On entry, ICM is an array of dimension 2. On exit, when JUMPM
 *         is  greater  than  zero,  this  array  contains  the  encoded
 *         constant  C  to jump from  X(n)  to  X(n+JUMPM) in the random
 *         sequence. ICM(0:1)  contains  respectively  the  16-lower and
 *         15-higher  bits  of this constant  C. When  JUMPM  is less or
 *         equal than zero, this array is not referenced.
 *
 * ---------------------------------------------------------------------
 */
   int                        j[2], k;

   if( JUMPM > 0 )
   {
      IAM[0] = MULT[0]; IAM[1] = MULT[1];   /* IAM   = MULT;          */
      ICM[0] = IADD[0]; ICM[1] = IADD[1];   /* ICM   = IADD;          */
      for( k = 1; k <= JUMPM-1; k++ )
      {
         ATL_lmul( IAM, MULT, j );          /* j     = IAM   * MULT;  */
         IAM[0] = j[0]; IAM[1] = j[1];      /* IAM   = j;             */
         ATL_lmul( ICM, MULT, j );          /* j     = ICM   * MULT;  */
         ATL_ladd( IADD, j, ICM );          /* ICM   = IADD  + j;     */
      }
      ATL_lmul( IRANN, IAM, j );            /* j     = IRANN * IAM;   */
      ATL_ladd( j, ICM, IRANM );            /* IRANM = j     + ICM;   */
   }
   else
   {                                        /* IRANM = IRANN          */
      IRANM[0] = IRANN[0]; IRANM[1] = IRANN[1];
   }
}


void ATL_srand(int iseed)
{
   int iadd[2], ia1[2], ic1[2], iran1[2], jseed[2], mult[2];

   mult [0] = 20077; mult[1] = 16838;
   iadd [0] = 12345; iadd [1] = 0;
   jseed[0] = iseed; jseed[1] = (iseed>>16);

   ATL_xjumpm( 1, mult, iadd, jseed, iran1, ia1, ic1 );
   ATL_setran( 0, iran1 ); ATL_setran( 1, ia1 ); ATL_setran( 2, ic1 );
}

int ATL_rand(void)
{
   int j[2];

   ATL_setran( 3, j );
   return(j[0] + ((j[1])<<16));
}
void printmat(char *mat, int M, int N, TYPE *A, int lda)
{
   int i, j;

#ifdef TCPLX
   lda *= 2;
#endif
   printf("\n%s = \n",mat);
   for (i=0; i != M; i++)
   {
#ifdef TREAL
      for (j=0; j != N; j++) printf("%f  ",A[i+j*lda]);
#else
      for (j=0; j != N; j++) printf("(%f,%f)  ",A[2*i+j*lda], A[1+2*i+j*lda]);
#endif
      printf("\n");
   }
}

void matgen(int M, int N, TYPE *A, int lda, int seed)
{
   int i, j;

#ifdef TCPLX
   M *= 2;
   lda *= 2;
#endif
   dumb_seed(seed);
   for (j=N; j; j--)
   {
      for (i=0; i != M; i++) A[i] = dumb_rand();
      A += lda;
   }
}

#define trusted_gemm(TA, TB, m, n, k, al, A, lda, B, ldb, be, C, ldc) \
   dgemmRCW(TA, TB, m, n, k, al, A, lda, B, ldb, be, C, ldc)

#define test_gemm(TA, TB, m, n, k, al, A, lda, B, ldb, be, C, ldc) \
   dgemm_asg(TA, TB, m, n, k, al, A, lda, B, ldb, be, C, ldc)

int mmcase(int TEST, int CACHESIZE, char TA, char TB, int M, int N, int K,
	   SCALAR alpha, TYPE *A, int lda, TYPE *B, int ldb, SCALAR beta,
	   TYPE *C, int ldc, TYPE *D, int ldd)
{
   char *pc;
#ifdef TREAL
   char *form="%4d   %c   %c %4d %4d %4d  %5.1f  %5.1f  %6.2f %5.1f %5.2f   %3s\n";
   #define MALPH alpha
   #define MBETA beta
#else
   #define MALPH *alpha, alpha[1]
   #define MBETA *beta, beta[1]
   char *form="%4d   %c   %c %4d %4d %4d  %5.1f %5.1f  %5.1f %5.1f  %6.2f %6.1f %4.2f   %3s\n";
#endif
   int ii, jj, i, j=0, PASSED, nerrs;
   double t0, t1, t2, t3, mflop;
   TYPE maxval, f1, ferr;
   static TYPE feps=0.0;
   static int itst=1;
   /*int *L2, nL2=(1.3*L2SIZE)/sizeof(int);*/
   enum ATLAS_TRANS TAc, TBc;
   double l2ret;

   if (!TEST) D = C;
   /*if (nL2) L2 = malloc(nL2*sizeof(int));*/
   l2ret = ATL_flushcache( CACHESIZE );
   if (TA == 'n' || TA == 'N')
   {
      matgen(M, K, A, lda, K*1112);
      TAc = AtlasNoTrans;
   }
   else
   {
      matgen(K, M, A, lda, K*1112);
      if (TA == 'c' || TA == 'C') TAc = AtlasConjTrans;
      else TAc = AtlasTrans;
   }
   if (TB == 'n' || TB == 'N')
   {
      matgen(K, N, B, ldb, N*2238);
      TBc = AtlasNoTrans;
   }
   else
   {
      matgen(N, K, B, ldb, N*2238);
      if (TB == 'c' || TB == 'C') TBc = AtlasConjTrans;
      else TBc = AtlasTrans;
   }
   matgen(M, N, C, ldc, M*N);

#ifdef DEBUG
   printmat("A0", M, K, A, lda);
   printmat("B0", K, N, B, ldb);
   printmat("C0", M, N, C, ldc);
#endif

   /*
     if (L2)
     {
     for (i=0; i != nL2; i++) L2[i] = 0.0;
     for (i=0; i != nL2; i++) j += L2[i];
     }*/

   /* invalidate L2 cache */
   l2ret = ATL_flushcache( -1 );

   t0 = time00();
   trusted_gemm(TAc, TBc, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   t1 = time00() - t0;
   if (t1 <= 0.0) mflop = t1 = 0.0;
   else   /* flop rates actually 8MNK+12MN & 2MNK + 2MN, resp */
      #ifdef TCPLX
         mflop = ( ((8.0*M)*N)*K ) / (t1*1000000.0);
      #else
         mflop = ( ((2.0*M)*N)*K ) / (t1*1000000.0);
      #endif
   printf(form, itst, TA, TB, M, N, K, MALPH, MBETA, t1, mflop, 1.0, "---");

#ifdef DEBUG
   printmat("C", M, N, C, ldc);
#endif

#ifndef TIMEONLY
   matgen(M, N, D, ldd, M*N);

   /* invalidate L2 cache */
   l2ret = ATL_flushcache( -1 );

   t0 = time00();
   test_gemm(TAc, TBc, M, N, K, alpha, A, lda, B, ldb, beta, D, ldd);

   t2 = time00() - t0;
//fprintf(stderr, "A=%f, B=%f, alpha=%f, good=%f, test=%f (%d)\n", *A, *B, alpha, *D,*C, C);
   if (t2 <= 0.0) t2 = mflop = 0.0;
   else
      #ifdef TCPLX
         mflop = ( ((8.0*M)*N)*K ) / (t2*1000000.0);
      #else
         mflop = ( ((2.0*M)*N)*K ) / (t2*1000000.0);
      #endif
#ifdef DEBUG
   printmat("D", M, N, D, ldd);
#endif
   if (TEST)
   {
      if (feps == 0.0)
      {
#if 0
         f1 = feps = 0.5;
         do
         {
            feps = f1;
            f1 *= 0.5;
            maxval = 1.0 + f1;
         }
         while (maxval != 1.0);
         printf("feps=%e\n",feps);
#else
         feps = EPS;
#endif
#ifdef DEBUG
         printf("feps=%e\n",feps);
#endif
      }
#ifdef TREAL
      ferr = 2.0 * (Mabs(alpha) * 2.0*K*feps + Mabs(beta) * feps) + feps;
#else
      f1 = Mabs(*alpha) + Mabs(alpha[1]);
      maxval = Mabs(*beta) + Mabs(beta[1]);
      ferr = 2.0 * (f1*8.0*K*feps + maxval*feps) + feps;
#endif
      PASSED = 1;
      maxval = 0.0;
      pc = "YES";
      nerrs = ii = jj = 0;
      for (j=0; j != N; j++)
      {
         for (i=0; i != M SHIFT; i++)
         {
            f1 = D[i] - C[i];
            if (f1 < 0.0) f1 = -f1;
            if (f1 > ferr)
            {
               nerrs++;
               PASSED = 0;
               pc = "NO!";
               if (f1 > maxval)
               {
                  maxval=f1;
                  ii = i+1;
                  jj = j+1;
               }
            }
         }
         D += ldd SHIFT;
         C += ldc SHIFT;
      }
      if (maxval != 0.0)
         fprintf(stderr, "ERROR: nerr=%d, i=%d, j=%d, maxval=%e\n", nerrs, ii,jj, maxval);
   }
   else pc = "---";
   if (t1 == t2) t3 = 1.0;
   else if (t2 != 0.0) t3 = t1/t2;
   else t3 = 0.0;
   printf(form, itst++, TA, TB, M, N, K, MALPH, MBETA, t2, mflop, t3, pc);
#else
   itst++;
   PASSED = 1;
#endif
   /*free(L2);*/
   l2ret = ATL_flushcache( 0 );
   return(PASSED);
}

int mmcase0(int MFLOP, int CACHESIZE, char TA, char TB, int M, int N, int K,
	    SCALAR alpha, int lda, int ldb, SCALAR beta, int ldc)
{
   char *pc;
#ifdef TREAL
   char *form="%4d   %c   %c %4d %4d %4d  %5.1f  %5.1f  %6.2f %5.1f %5.2f   %3s\n";
   #define MALPH alpha
   #define MBETA beta
   TYPE betinv, bet=beta;
#else
   #define MALPH *alpha, alpha[1]
   #define MBETA *beta, beta[1]
   char *form="%4d   %c   %c %4d %4d %4d  %5.1f %5.1f  %5.1f %5.1f  %6.2f %6.1f %4.2f   %3s\n";
   TYPE betinv[2], *bet=beta;
#endif
   int nreps, incA, incB, incC, inc, nmat, k;
   TYPE *c, *C, *a, *A, *b, *B, *st;
   int ii, jj, i, j=0, PASSED, nerrs;
   double t0, t1, t2, t3, mflop, mf, mops;
   TYPE maxval, f1, ferr;
   static TYPE feps=0.0;
   static int itst=1;
   enum ATLAS_TRANS TAc, TBc;
   void *vp;

   #ifdef TCPLX
      if (*beta == 0.0 && beta[1] == 0.0) betinv[0] = betinv[1] = 0.0;
      else if (beta[1] == 0.0) { betinv[0] = 1 / *beta;  betinv[1] = 0.0; }
      else
      {
         t0 = *beta;
         t1 = beta[1];
         if (Mabs(t1) <= Mabs(t0))
         {
            t2 = t1 / t0;
            betinv[0] = t0 = 1.0 / (t0 + t1*t2);
            betinv[1] = -t0 * t2;
         }
         else
         {
            t2 = t0 / t1;
            betinv[1] = t0 = -1.0 / (t1 + t0*t2);
            betinv[0] = -t2 * t0;
         }
      }
      mops = ( ((8.0*M)*N)*K ) / 1000000.0;
   #else
      if (beta != 0.0) betinv = 1.0 / beta;
      else betinv = beta;
      mops = ( ((2.0*M)*N)*K ) / 1000000.0;
   #endif
   nreps = MFLOP / mops;
   if (nreps < 1) nreps = 1;
   if (TA == 'n' || TA == 'N')
   {
      TAc = AtlasNoTrans;
      incA = lda * K;
   }
   else
   {
      if (TA == 'c' || TA == 'C') TAc = AtlasConjTrans;
      else TAc = AtlasTrans;
      incA = lda * M;
   }
   if (TB == 'n' || TB == 'N')
   {
      incB = ldb * N;
      TBc = AtlasNoTrans;
   }
   else
   {
      incB = ldb * K;
      if (TB == 'c' || TB == 'C') TBc = AtlasConjTrans;
      else TBc = AtlasTrans;
   }
   incC = ldc*N;
   inc = incA + incB + incC;
   i = M*K + K*N + M*N;  /* amount of inc actually referenced */
   /* This is a hack; change to use of flushcache instead. */
   nmat = ((CACHESIZE/ATL_sizeof) + i)/i;
   vp = malloc(ATL_MulBySize(nmat*inc)+ATL_Cachelen);
   ATL_assert(vp);
   C = c = ATL_AlignPtr(vp);
   a = A = C + incC;
   b = B = A + incA;
   st = C + nmat*inc;
   matgen(inc, nmat, C, inc, M*N);

#ifdef DEBUG
   printmat("A0", M, K, A, lda);
   printmat("B0", K, N, B, ldb);
   printmat("C0", M, N, C, ldc);
#endif

   t0 = time00();
   for (k=nreps; k; k--)
   {
      trusted_gemm(TAc, TBc, M, N, K, alpha, a, lda, b, ldb, bet, c, ldc);
      c += inc; a += inc; b += inc;
      if (c == st)
      {
         c = C; a = A; b = B;
         if (bet == beta) bet = betinv;
         else bet = beta;
      }
   }
   t1 = time00() - t0;
   t1 /= nreps;
   if (t1 <= 0.0) mflop = t1 = 0.0;
   else   /* flop rates actually 8MNK+12MN & 2MNK + 2MN, resp */
      mflop = mops / t1;
   printf(form, itst, TA, TB, M, N, K, MALPH, MBETA, t1, mflop, 1.0, "---");

#ifdef DEBUG
   printmat("C", M, N, C, ldc);
#endif

   matgen(inc, nmat, C, inc, M*N);
   t0 = time00();
   for (k=nreps; k; k--)
   {
      test_gemm(TAc, TBc, M, N, K, alpha, a, lda, b, ldb, bet, c, ldc);
      c += inc; a += inc; b += inc;
      if (c == st)
      {
         c = C; a = A; b = B;
         if (bet == beta) bet = betinv;
         else bet = beta;
      }
   }

   t2 = time00() - t0;
   t2 /= nreps;
   if (t2 <= 0.0) t2 = mflop = 0.0;
   else mflop = mops / t2;

   pc = "---";
   if (t1 == t2) t3 = 1.0;
   else if (t2 != 0.0) t3 = t1/t2;
   else t3 = 0.0;
   printf(form, itst++, TA, TB, M, N, K, MALPH, MBETA, t2, mflop, t3, pc);
   free(vp);
   return(1);
}

void PrintUsage(char *nam)
{
   fprintf(stderr, "USAGE: %s -# <nreps> -Atrans <ntrans> n/t/c -Btrans <ntrans> n/t/c -M <m1> <mN> <minc> -N <n1> <nN> <ninc> <k1> <kN> <kinc> -n <n> -m <m> -k <k> -a <nalphas> <alpha1> ... <alphaN> -b <nbetas> <beta1> ... <betaN> -Test <0/1> -F <mflops> -C <cachesize>\n", nam);
   exit(-1);
}

void GetFlags(int nargs, char *args[], int *nrep, int *TEST,
              int *nta, enum ATLAS_TRANS **TransA,
              int *ntb, enum ATLAS_TRANS **TransB,
              int *M0, int *MN, int *Minc, int *N0, int *NN, int *Ninc, 
              int *K0, int *KN, int *Kinc,
              int *nalphas, TYPE **alphas, int *nbetas, TYPE **betas,
              int *LDA_IS_M, int *MFLOP, int *CACHESIZE)

{
   char ch;
   int i, j;
/*
 * Set up defaults
 */
   *nrep = 1;
   *TEST = 1;
   *M0 = *N0 = *K0 = -1;
   *nta = *ntb = *nalphas = *nbetas = -1;
   *MFLOP = *LDA_IS_M = 0;
#ifdef L2SIZE
   *CACHESIZE = L2SIZE;               /* Size of largest cache to flush */
#else
   *CACHESIZE = 4*1024*1024;
#endif

   for (i=1; i < nargs; i++)
   {
      switch(args[i][1])
      {
      case 'F':
         *MFLOP = atoi(args[++i]);
         break;
      case 'C':
            if( args[i+1] == NULL ) PrintUsage( args[0] );
	    *CACHESIZE = 1024*atoi(args[++i]);
            break;
      case '#':
         *nrep = atoi(args[++i]);
         break;
      case 'A':
         *nta   = atoi(args[++i]);
         *TransA = malloc(*nta * sizeof(int));
         assert(*TransA);
         for (j=0; j != *nta; j++)
         {
            ch = *args[++i];
            if (ch == 'n' || ch == 'N') (*TransA)[j] = AtlasNoTrans;
            else if (ch == 't' || ch == 'T') (*TransA)[j] = AtlasTrans;
            else if (ch == 'c' || ch == 'C') (*TransA)[j] = AtlasConjTrans;
            else PrintUsage(args[0]);
         }
         break;
      case 'B':
         *ntb   = atoi(args[++i]);
         *TransB = malloc(*ntb * sizeof(int));
         assert(*TransB);
         for (j=0; j != *ntb; j++)
         {
            ch = *args[++i];
            if (ch == 'n' || ch == 'N') (*TransB)[j] = AtlasNoTrans;
            else if (ch == 't' || ch == 'T') (*TransB)[j] = AtlasTrans;
            else if (ch == 'c' || ch == 'C') (*TransB)[j] = AtlasConjTrans;
            else PrintUsage(args[0]);
         }
         break;
      case 'M':
         *M0 = atoi(args[++i]);
         *MN = atoi(args[++i]);
         *Minc = atoi(args[++i]);
         break;
      case 'N':
         *N0 = atoi(args[++i]);
         *NN = atoi(args[++i]);
         *Ninc = atoi(args[++i]);
         break;
      case 'K':
         *K0 = atoi(args[++i]);
         *KN = atoi(args[++i]);
         *Kinc = atoi(args[++i]);
         break;
      case 'T':
         *TEST = atoi(args[++i]);
         break;
      case 'm':
         *M0 = *MN = *Minc = atoi(args[++i]);
         break;
      case 'n':
         *N0 = *NN = *Ninc = atoi(args[++i]);
         break;
      case 'k':
         *K0 = *KN = *Kinc = atoi(args[++i]);
         break;
      case 'a':
         *nalphas = atoi(args[++i]);
         *alphas = malloc(ATL_MulBySize(*nalphas));
         assert(*alphas);
         for (j=0; j < *nalphas SHIFT; j++) (*alphas)[j] = atof(args[++i]);
         break;
      case 'b':
         *nbetas  = atoi(args[++i]);
         *betas  = malloc(ATL_MulBySize(*nbetas ));
         assert(*betas );
         for (j=0; j < *nbetas SHIFT; j++) (*betas)[j] = atof(args[++i]);
         break;
      case 'd':
         *LDA_IS_M  = atoi(args[++i]);
         break;
      default:
         PrintUsage(args[0]);
         break;
      }
   }
/*
 * Finish setting up defaults if the user has not selected
 */
   if (*N0 == -1)
   {
      *N0 = 100;
      *NN = 1000;
      *Ninc = 100;
   }
   if (*nta == -1)
   {
      *nta = 1;
      *TransA = malloc(sizeof(int));
      assert(*TransA);
      **TransA = AtlasTrans;
   }
   if (*ntb == -1)
   {
      *ntb = 1;
      *TransB = malloc(sizeof(int));
      assert(*TransB);
      **TransB = AtlasNoTrans;
   }
   if (*nalphas == -1)
   {
      *nalphas = 1;
      *alphas = malloc(ATL_MulBySize(1));
      assert(*alphas);
      #ifdef TREAL
         **alphas = 1.0;
      #else
         **alphas = 1.0;
         (*alphas)[1] = 0.0;
      #endif
   }
   if (*nbetas  == -1)
   {
      *nbetas  = 1;
      *betas  = malloc(ATL_MulBySize(1));
      assert(*betas );
      #ifdef TREAL
         **betas  = 1.0;
      #else
         **betas  = 1.0;
         (*betas)[1] = 0.0;
      #endif
   }
}
___main(){}
__main(){}
MAIN__(){}
_MAIN_(){}
main(int nargs, char *args[])
/*
 *  tst <tst> <# TA> <TA's> <# TB's> <TB's> <M0> <MN> <incM> <N0> <NN> <incN>
 *      <K0> <KN> <incK> <# alphas> <alphas> <# betas> <betas>
 *
 */
{
   int M0, MN, incM, N0, NN, incN, K0, KN, incK, lda, ldb, ldc, MFLOP;
   int i, k, m, n, im, in, ik, ita, itb, ia, ib, nTA, nTB, nalph, nbeta;
   int ir, nrep;
   int itst=0, ipass=0, TEST, LDA_IS_M, MSAME=0, KSAME=0;
   TYPE *alph, *beta, *A, *B, *C, *D=NULL;
   #ifdef TREAL
      TYPE bet1 = 1.0, alp1 = -1.0;
   #else
      TYPE bet1[2] = {1.0, 0.0}, alp1[2] = {-1.0, 0.0};
   #endif
   char TA, TB;
   enum ATLAS_TRANS *TransA, *TransB, TAc, TBc;
   int CACHESIZE;

   GetFlags(nargs, args, &nrep, &TEST, &nTA, &TransA, &nTB, &TransB, 
            &M0, &MN, &incM, &N0, &NN, &incN, &K0, &KN, &incK,
            &nalph, &alph, &nbeta, &beta, &LDA_IS_M, &MFLOP, &CACHESIZE);

   if (M0 == -1)
   {
      MSAME = 1;
      M0 = MN = incM = NN;
   }
   if (K0 == -1)
   {
      KSAME = 1;
      K0 = KN = incK = NN;
   }

   if (!MFLOP)
   {
      A = malloc(MN*KN*ATL_sizeof);
      B = malloc(NN*KN*ATL_sizeof);
      C = malloc(MN*NN*ATL_sizeof);
      if (TEST) D = malloc(MN*NN*ATL_sizeof);
      else D = NULL;
      if (!A || !B || !C || (TEST && !D))
      {
         fprintf(stderr, "Not enough memory to run tests!!\n");
         exit(-1);
      }
   }
/*
 * Page the code in from disk, so first timing doesn't blow
 */
#if 0
   if (MFLOP)
   {
      mmcase0(10, 1, 'n', 'n', 100, 100, 100, alp1, 100, 100, bet1, 100);
      mmcase0(10, 1, 'n', 't', 100, 100, 100, alp1, 100, 100, bet1, 100);
      mmcase0(10, 1, 't', 'n', 100, 100, 100, alp1, 100, 100, bet1, 100);
      mmcase0(10, 1, 't', 't', 100, 100, 100, alp1, 100, 100, bet1, 100);
   }
   else
   {
      m = Mmin(100, MN);
      k = Mmin(100, KN);
      n = Mmin(100, NN);
      matgen(m, k, A, m, m*k);
      matgen(k, n, B, k, n*k);
      matgen(m, n, C, m, m*n);
      TA = TB = 'N';
      TAc = TBc = AtlasNoTrans;
      trusted_gemm(TAc, TBc, m, n, k, alp1, A, m, B, k, bet1, C, m);
      test_gemm(TAc, TBc, m, n, k, alp1, A, m, B, k, bet1, C, m);
   }
#endif

#ifdef TREAL
   printf("\nTEST  TA  TB    M    N    K  alpha   beta    Time  Mflop  SpUp  PASS\n");
   printf("====  ==  ==  ===  ===  ===  =====  =====  ======  =====  ====  ====\n\n");
#else
   printf("\nTEST  TA  TB    M    N    K        alpha         beta    Time  Mflop  SpUp  PASS\n");
   printf("====  ==  ==  ===  ===  ===  ===== =====  ===== =====  ======  =====  ====  ====\n\n");
#endif
   for (im=M0; im <= MN; im += incM)
   {
      for (n=N0; n <= NN; n += incN)
      {
         if (MSAME) m = n;
         else m = im;
         for (ik=K0; ik <= KN; ik += incK)
         {
            if (KSAME) k = n;
            else k = ik;
            for (ita=0; ita != nTA; ita++)
            {
               if (TransA[ita] == AtlasNoTrans) TA = 'N';
               else if (TransA[ita] == AtlasTrans) TA = 'T';
               else if (TransA[ita] == AtlasConjTrans) TA = 'C';

               for (itb=0; itb != nTB; itb++)
               {
                  if (TransB[itb] == AtlasNoTrans) TB = 'N';
                  else if (TransB[itb] == AtlasTrans) TB = 'T';
                  else if (TransB[itb] == AtlasConjTrans) TB = 'C';
                  if (LDA_IS_M)
                  {
                     if (TA == 'n' || TA == 'N') lda = m;
                     else lda = k;
                     if (TB == 'n' || TB == 'N') ldb = k;
                     else ldb = n;
                     ldc = m;
                  }
                  else
                  {
                     if (TA == 'n' || TA == 'N') lda = MN;
                     else lda = KN;
                     if (TB == 'n' || TB == 'N') ldb = KN;
                     else ldb = NN;
                     ldc = MN;
                  }
                  for (ia=0; ia != nalph; ia++)
                  {
                     for (ib=0; ib != nbeta; ib++)
                     {
                        for (ir=0; ir < nrep; ir++)
                        {
                           itst++;
                           if (MFLOP)
                           {
                              ipass++;
#ifdef TREAL
                                 mmcase0(MFLOP, CACHESIZE, TA, TB, m, n, k,
				         alph[ia], lda, ldb, beta[ib], ldc);
#else
                                 mmcase0(MFLOP, CACHESIZE, TA, TB, m, n, k,
				         alph+(ia SHIFT), lda, ldb,
				         beta+(ib SHIFT), ldc);
#endif
                           }
                           else
                           {
#ifdef TREAL
                                 ipass += mmcase(TEST, CACHESIZE, TA, TB, m,
					         n, k, alph[ia], A, lda, B, ldb,
                                                 beta[ib], C, ldc, D,ldc);
#else
                                 ipass += mmcase(TEST, CACHESIZE, TA, TB, m,
					         n, k, alph+(ia SHIFT), A,
					         lda, B, ldb, beta+(ib SHIFT),
					         C, ldc, D,ldc);
#endif
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   if (TEST && !MFLOP)
      printf("\nNTEST=%d, NUMBER PASSED=%d, NUMBER FAILURES=%d\n",
             itst, ipass, itst-ipass);
   else printf("\nDone with %d timing runs\n",itst);
   free(TransA);
   free(TransB);
   free(alph);
   free(beta);
   if (!MFLOP)
   {
      free(A);
      free(B);
      free(C);
      if (D) free(D);
   }
   exit(0);
}
