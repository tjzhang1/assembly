#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

unsigned long myGap(void *p0, void *p1)
{
   size_t a0=(size_t)p0, a1=(size_t)p1;
   return((a0 <= a1) ? a1-a0 : a0-a1);
}

typedef struct items item_t;
struct items
{
   char name[32];       /* name of item */
   item_t *next;        /* ptr to next item in queue */
   float cost, weight;  /* mul by num for weight of all items */
   int i0, i1;          /* usually low & hi damage, but varies by item */
   unsigned int num;    /* number of this item being carried */
   unsigned int dsc;    /* if non-zero, extended description */
   unsigned int flg;    /* bitvec used to store info about item */
};

item_t *findMini0(item_t *p)
{
   item_t *minp=p;
   if (p)
   {
      int minI=p->i0;
      p = p->next;
      while (p)
      {
         int i = p->i0;
         if (minI > i)
         {
            minI = i;
            minp = p;
         }
         p = p->next;
      }
   }
   return(minp);
}

item_t *findMini1(item_t *p)
{
   item_t *minp=p;
   if (p)
   {
      int minI=p->i1;
      p = p->next;
      while (p)
      {
         int i = p->i1;
         if (minI > i)
         {
            minI = i;
            minp = p;
         }
         p = p->next;
      }
   }
   return(minp);
}

void chkItem(char wht, item_t *bp, item_t *prev, item_t *good, item_t *test)
{
   if (good != test)
   {
      fprintf(stderr, "ERROR %u: you get minp=%p, correct=%p\n", __LINE__,
              test, good);
      if (wht == '0')
         fprintf(stderr, "   YOUR IVAL=%d (correct=%d)\n", 
                 test->i0, good->i0);
      else
         fprintf(stderr, "   YOUR IVAL=%d (correct=%d)\n", 
                 test->i1, good->i1);
      exit(10);
   }
   if (prev)
   {
      if (prev->next != good)
      {
         fprintf(stderr, "PREV BAD prev=%p, prev->next=%p, correct next=%p\n",
                 prev, prev->next, good);
         exit(11);
      }
   }
   else
   {
      if (good != bp)
         fprintf(stderr, "PREV=NULL, but min=%p and base=%p!\n", good, bp);
      assert(good == bp);
   }
   fprintf(stderr, "PASS BASIC ITEM i%c SEARCH\n", wht);
}

typedef struct word word_t;
struct word
{
   char *word;
   word_t *next;
   int len;        // len of word, not including str term
};

word_t *findMinLen(word_t *p)
{
   word_t *minp=p;
   if (p)
   {
      int minLen=p->len;
      p = p->next;
      while (p)
      {
         int len = p->len;
         if (minLen > len)
         {
            minLen = len;
            minp = p;
         }
         p = p->next;
      }
   }
   return(minp);
}

int main(int nargs, char **args)
{
   unsigned int i, N=17, nerr=0;
   unsigned long gapG, gapT, gapI;
   void *vp;
   item_t *ip, *ib=NULL, *imin;
   word_t *wp, *wb=NULL, *wmin;
   unsigned long addrGap(void *a0, void *a1);
   int cntNodes(void *bp, unsigned long gapNxt);
   void *findIntMinNode(void *bp, void *minPrev, 
                        unsigned long nextgap, unsigned long igap);

   if (nargs > 1)
   {
      N = atoi(args[1]);
      assert(N < 5000 && N > 1);
   } 
   gapG = myGap(wp, &(wp->next));
   gapT = addrGap(wp, &(wp->next));
   if (gapG != gapT)
   {
      fprintf(stderr, "addGap ERROR: exp=%lu, got=%lu\n", gapG, gapT);
      nerr++;
   }
   gapT = addrGap(&(wp->next), wp);
   if (gapG != gapT)
   {
      fprintf(stderr, "addGap REVERSE ERROR: exp=%lu, got=%lu\n", gapG, gapT);
      nerr++;
   }
   gapG = myGap(ip, &(ip->next));
   gapT = addrGap(ip, &(ip->next));
   if (gapG != gapT)
   {
      fprintf(stderr, "addGap ITEM ERROR: exp=%lu, got=%lu\n", gapG, gapT);
      nerr++;
   }
   if (nerr)       // Don't test cnt until addrGap works
      return(nerr);
/*
 * Test if cntNodes handles NULL case
 */
   i = cntNodes(NULL, -1);
   if (i != 0)
   {
      fprintf(stderr, "you cannot count to 0!\n");
      return(1);
   }
/*
 * Allocate N-length lists of words and items, with random int vals
 */
   srand(N);
   for (i=0; i < N; i++)
   {
      ip = malloc(sizeof(item_t));
      assert(ip);
      wp = malloc(sizeof(word_t));
      assert(wp);
      wp->next = wb;
      wb = wp;
      ip->next = ib;
      ib = ip;
      wp->len = rand();
      ip->i0 = rand();
      ip->i1 = rand();
   }
/*
 * See if you can count to 1 only
 */
   gapT = addrGap(ip, &(ip->next));
   gapG = addrGap(wp, &(wp->next));
   vp = ib->next;
   ib->next = NULL;
   assert(cntNodes(ib, gapT) == 1);
   ib->next = vp;

   vp = wb->next;
   wb->next = NULL;
   assert(cntNodes(wb, gapG) == 1);
   wb->next = vp;
   fprintf(stderr, "You can count to 1 (birds can count to 5) . . .\n");

   i = cntNodes(wb, gapG);
   if (i != N)
      fprintf(stderr, "ERROR %u, line %u: exp=%u, got=%u\n", ++nerr, __LINE__,
              N, i);

   i = cntNodes(ib, gapT);
   if (i != N)
      fprintf(stderr, "ERROR %u, line %u: exp=%u, got=%u\n", ++nerr, __LINE__,
              N, i);
   if (nerr)       // Don't test findMin until cntNodes works
      return(nerr);


   fprintf(stderr, "TRYING NULL SEARCH\n");
   wp = (wb->next) ? wb->next: wb;          // init wp to non-NULL value
   vp = findIntMinNode(NULL, &wp, -1, -1);  // see if NULL case handled
   assert(vp == NULL);
   assert(wp == NULL);
   fprintf(stderr, "PASS   NULL SEARCH!\n");

   fprintf(stderr, "\nTRYING NO-PREV NULL SEARCH!\n");
   vp = findIntMinNode(NULL, NULL, -1, -1);  // see if NULL case handled
   fprintf(stderr, "PASS   NO-PREV NULL SEARCH!\n\n");

   gapI = addrGap(wb, &(wb->len));
   vp = findIntMinNode(wb, &wp, gapG, gapI);
   wmin = findMinLen(wb);
   if (wmin != vp)
   {
      fprintf(stderr, "ERROR %u: you get minp=%p, correct=%p\n", __LINE__,
              vp, wmin);
      fprintf(stderr, "   YOUR IVAL=%d (correct=%d)\n", 
              ((word_t*)vp)->len, wmin->len);
      return(1);
   }
   if (wp)
      assert(wp->next == wmin);
   else
      assert(wmin == wb);
   fprintf(stderr, "PASS BASIC WORD SEARCH\n");
   while (wb)
   {
      wp = wb->next;
      free(wb);
      wb = wp;
   }

   gapI = addrGap(ib, &(ib->i0));
   vp = findIntMinNode(ib, &ip, gapT, gapI);
   imin = findMini0(ib);
//   fprintf(stderr, "   bp=%p, gapI=%lu, min=%p,%p\n", ib, gapI, imin, vp);
   chkItem('0', ib, ip, imin, vp);
   fprintf(stderr, "TRYING WITH NO PREV:\n");
   vp = findIntMinNode(ib, NULL, gapT, gapI);
   assert(vp == imin);
   fprintf(stderr, "PASS   WITH NO PREV:\n");

/*
 * Force min to be first node
 */
   ib->i0 = imin->i0 - 1;
   assert(ib->i0 < imin->i0);  /* assert no underflow */
   vp = findIntMinNode(ib, &ip, gapT, gapI);
   assert(ip == NULL);
   assert(vp == ib);

   vp = findIntMinNode(ib, &ip, gapT, addrGap(ib, &(ib->i1)));
   imin = findMini1(ib);
//   fprintf(stderr, "   imin=%p,%p, prev=%p, bp=%p\n", imin, vp, ip, ib);
   chkItem('1', ib, ip, imin, vp);

   while (ib)
   {
      ip = ib->next;
      free(ib);
      ib = ip;
   }

   if (!nerr)
      fprintf(stderr, "PASSED ALL TESTS!\n");
   else
      fprintf(stderr, "FAILED %u TESTS!\n", nerr);
   return(0);
}
