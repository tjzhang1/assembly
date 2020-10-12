#include <stdio.h>
#include <string.h>
#include <assert.h>
int main(void)
{
   int iv = 5, ires;
   char alph[28];
   #ifndef IGNORE_STR
      char *exp = "abcdefghijklmnopqrstuvwxyz";
   #endif
   int istr(int i1, int i2, int i3, int *i4, char *alp);

   ires = istr(4, 3, 7, &iv, alph);
   if (ires != 7)
   {
      fprintf(stderr, "ires should be=%d; got=%d\n", 7, ires);
      assert(0);
   }
   #ifndef IGNORE_IV
      if (iv != 5*5+7)
      {
         fprintf(stderr, "IV: expected=%d, got=%d\n", 5*5+7, iv);
         assert(0);
      }
   #endif
   #ifndef IGNORE_STR
      if (alph[26] != '\0')
         fprintf(stderr, "Don't see string term at end of string!\n");
      if (strcmp(alph, exp))
      {
         fprintf(stderr, "STR: expected='%s'\n", exp);
         fprintf(stderr, "STR:      got='%s'\n", alph);
         assert(0);
      }
   #endif

   ires = istr(11, 17, 3,  &iv, alph);
   if (ires != 28)
   {
      fprintf(stderr, "ires should be=%d; got=%d\n", 28, ires);
      assert(0);
   }
   #ifndef IGNORE_IV
      if (iv != (5*5+7)*5+3)
      {
         fprintf(stderr, "IV: expected=%d, got=%d\n", (5*3+7)*5+3, iv);
         assert(0);
      }
   #endif
   #ifndef IGNORE_STR
      if (strcmp(alph, exp))
      {
         fprintf(stderr, "STR: expected='%s'\n", exp);
         fprintf(stderr, "STR:      got='%s'\n", alph);
         assert(0);
      }
   #endif
   #if defined(IGNORE_STR) && defined(IGNORE_IV)
      printf("Throw a party, you can add two integers.  How about pointers?\n");
   #elif defined(IGNORE_STR)
      printf("OK, integers work, but what about strings?\n");
   #else
      printf("Flawless victory!\n");
   #endif
   return(0);
}
