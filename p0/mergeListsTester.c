#include <stdlib.h>
#include <stdio.h>
#include "word.h"
#include <string.h>
#include "mrgSort.h"

int main(void) {
   int i;
   word_t *a, *b, *curr;

   //init a
   a = malloc(sizeof(word_t));
   a->word = "hello";
   a->len = 1;
   curr = a;
   for(i=12; i>2; i-=2) {
      curr->next = malloc(sizeof(word_t));
      curr = curr->next;
      curr->word = "hello";
      curr->len = i;
   }
   curr->next = NULL;

   //init b
   b = malloc(sizeof(word_t));
   b->word = "hello";
   b->len = 25;
   curr = b;
   for(i=7; i<9; i++) {
      curr->next = malloc(sizeof(word_t));
      curr = curr->next;
      curr->word = "hello";
      curr->len = i;
   }
   curr->next = NULL;

   //print a
   curr = a;
   printf("a: [");
   while(curr) {
      printf("%d, ", curr->len);
      curr = curr->next;
   }
   printf("]\n");

   //print b
   curr = b;
   printf("b: [");
   while(curr) {
      printf("%d, ", curr->len);
      curr = curr->next;
   }
   printf("]\n");

   a = mergeLists(6, 3, a, b);
   //print a
   curr = a;
   printf("a: [");
   while(curr) {
      printf("%d, ", curr->len);
      curr = curr->next;
   }
   printf("]\n");

   a = mrgSortFast(a);
   curr = a;
   printf("a: [");
   while(curr) {
      printf("%d, ", curr->len);
      curr = curr->next;
   }
   printf("]\n");
   return 0;

}
