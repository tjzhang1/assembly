#ifndef WORD_H    /* guard multiple includes */
   #define WORD_H
typedef struct word word_t;
struct word
{
   word_t *next;
   #ifndef SMALL_WORD
   char word[128];
   #endif
   int len;        // len of word, not including str term
};
#endif
