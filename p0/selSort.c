#include <stdlib.h>
#include "word.h"

word_t *findMinPrev(word_t *lb)
/*
 * RETURNS: address of node PREVIOUS to node holding minimum ->len value.
 *          If list empty or base has minimum, return NULL.
 *          If two nodes have the same minimum ->len, return 1st.
 */
{
   if (lb) /* don't seg fault for lb==NULL */
   {
      if (lb->next) /* have at least 2 nodes, need to search */
      {
         register int imin=lb->len;
         word_t *p, *prev=lb, *pret=NULL;

         for (p=lb->next; p; p = p->next)
         {
            register int val=p->len;
            if (val < imin)
            {
               pret = prev;
               imin = val;
            }
            prev = p;
         }
         return(pret);
      }  /* fall thru to NULL return if lb only node so has min */
   }
   return(NULL);
}

word_t *selSort(word_t *ub)
/*
 * On entry ub is an unsorted list, which will be destroyed to produce sort.
 * RETURNS: base ptr of list sorted by ->len.
 */
{
   if (ub)
   {
      if (ub->next) /* need more than one node, or already sorted */
      {
         word_t *sb, *p, *last;
/*       Initialize sb outside loop to avoid extra ifs inside */
         p = findMinPrev(ub);
         if (p)
         {
            sb = p->next;       /* sorted base is minimum node */
            p->next = sb->next; /* remove minnode frm unsorted list */
         }
         else  /* prev==NULL means minnode is unsorted base node */
         {
            sb = ub;            /* sorted base is minimum node */
            ub = ub->next;      /* remove minnode frm unsorted list */
         }
         last = sb;             /* last node is 1st for one-length list */
         while(ub)
         {
            p = findMinPrev(ub);
            if (p)                      /* minnode is internal to ub */
            {
               last->next = p->next;    /* link sorted list to min node */
               p->next = p->next->next; /* remove min node from old list */
            }
            else                        /* minnode is ub */
            {
               last->next = ub;    /* link sorted list to min node */
               ub = ub->next;      /* remove min node from old list */
            }
            last = last->next;     /* move last to newly added sorted node */
         }
         last->next = NULL;       /* make sorted list NULL terminated */
         return(sb);
      }
      /* 1-length list will fall thru to returning unsorted list */
   }
   return(ub);  /* Zero or one length list already sorted */
}
