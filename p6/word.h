#ifndef WORD_H    /* guard multiple includes */
   #define WORD_H
typedef struct word word_t;
struct word
{
   char *word;
   word_t *next;
   int len;        // len of word, not including str term
};
#endif
