#ifndef ATLAS_DEF_H
   #define ATLAS_DEF_H

#include <stdio.h>
#include <stdlib.h>
#define time00 ATL_cputime
#ifndef TYPE
   #define TREAL
   #ifdef DREAL
      #define EPS 1.0e-15
      #define ATL_depsilon EPS
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
      #define ATL_sizeof 8
      #define ATL_MulBySize(i_) ((i_)<<3)
      #define ATL_DivBySize(i_) ((i_)>>3)
   #else
      #define EPS 5.0e-7 
      #define TYPE float
      #define PRE s
      #define UPR s
      #define PREU S
      #define PATL ATL_s
      #define PATU ATLU_s
      #define UATL ATLU_s
      #define CBLA cblas_s
      #define PATLU ATL_s
      #define ATL_rone   1.0f
      #define ATL_rnone -1.0f
      #define ATL_rzero   0.0f
      #define ATL_typify(m_) Mjoin(m_,f)
      #define SHIFT
      #define SCALAR TYPE
      #define SADD &
      #define SVAL
      #define SVVAL *
      #define ATL_sizeof 4
      #define ATL_MulBySize(i_) ((i_)<<2)
      #define ATL_DivBySize(i_) ((i_)>>2)
   #endif
   #define SCALAR_IS_ONE(M_scalar) ((M_scalar) == ATL_rone)
   #define SCALAR_IS_NONE(M_scalar) ((M_scalar) == ATL_rnone)
   #define SCALAR_IS_ZERO(M_scalar) ((M_scalar) == ATL_rzero)
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

double ATL_cputime(void);
#endif
