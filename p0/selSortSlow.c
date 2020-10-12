#include <stdlib.h>
#include "word.h"

word_t *findMin(word_t *lb)
/*
 * RETURNS: address of node with minimum ->len in lb
 */
{
   if (lb) /* don't seg fault for lb==NULL */
   {
      register int imin=lb->len;
      word_t *minp=lb, *p;

      for (p=lb->next; p; p = p->next)
      {
         register int val=p->len;
         if (val < imin)
         {
            minp = p;
            imin = val;
         }
      }
      return(minp);
   }
   return(NULL);
}

word_t *findPrev(word_t *lb, word_t *targ)
/*
 * RETURNS: address of node previous to targ in list starting at lb
 *          Will return NULL if lb or targ is NULL, or if targ not found.
 */
{
   word_t *p, *prev=NULL;
   for (p=lb; p; p = p->next)
   {
      if (p == targ)
         return(prev);
      prev = p;
   }
   return(prev);
}

word_t *selSortSlow(word_t *ub)
/*
 * ub is the base of unsorted list.  Other two params are ignored so this
 * can be a classic selection sort, which is O(N^2).
 *
 * On entry ub is an unsorted list, which will be destroyed to produce sort.
 * RETURNS: base ptr of list sorted by ->len.
 */
{
   word_t *sb, *last, *min, *prev;
   if (!ub)            /* empty list already sorted */
      return(NULL);    /* return NULL */
   if (!ub->next)      /* 1-entry list already sorted */
      return(ub);      /* return original list */
/*       Initialize sb outside loop to avoid extra ifs inside */
   min = findMin(ub);
   prev = findPrev(ub, min);
   if (prev)
   {
      sb = min;              /* sorted base is minimum node */
      prev->next = sb->next; /* remove minnode frm unsorted list */
   }
   else  /* prev==NULL means minnode is unsorted base node */
   {
      sb = ub;            /* sorted base is minimum node */
      ub = ub->next;      /* remove minnode frm unsorted list */
   }
   last = sb;             /* last node is 1st for one-length list */
   while(ub)
   {
      min = findMin(ub);
      prev = findPrev(ub, min);
      if (prev)                   /* minnode is internal to ub */
      {
         last->next = min;        /* link sorted list to min node */
         prev->next = min->next;  /* remove min node from old list */
      }
      else                        /* minnode is ub */
      {
         last->next = ub;         /* link sorted list to min node */
         ub = ub->next;           /* remove min node from old list */
      }
      last = last->next;          /* move last to newly added sorted node */
   }
   last->next = NULL;             /* make sorted list NULL terminated */
   return(sb);
}
