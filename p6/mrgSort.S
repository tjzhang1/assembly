
/*
r0                      r0               r1                 r2
word_t *mrgSort(word_t *bp, unsigned int gNxt, unsigned int gI)
{
   if(bp) //ensure there's at least 1 element
   {
      register int N=1, prevLen=bp->len, sorted=1;
      word_t *curr, *UB=NULL, *ML=NULL;

      //find length of list and determine whether it's sorted already
      for(curr=bp->next; curr; curr=curr->next)
      {
         register int currLen = curr->len;
         sorted &= (prevLen <= currLen); //remains 1 if sorted, 0 if not
         /update prevLen before incrementing curr
         prevLen = currLen;
         N++;
      }
      //if found to be sorted, return
      if(sorted)
         return(bp);
      //otherwise, sort
      return mrgSortR(N,bp,gNxt,gI,&ML,&UB);
   }
   return bp;
}
*/

.text
.thumb
.syntax unified
.local LOOP,DONE,QRET1,QRET2
.global mrgSort

#define bp r0
#define gNxt r1
#define gI r2
#define N r3
#define pML r4
#define pUB r5
#define currLen r6
#define curr r7

mrgSort:
   cmp bp,#0           //if no elements
   bEQ QRET1           //   return bp
   push {r4-r7,LR}     //save caller's variables
   mov r4,r8           //wrapper only uses r4-r7, but mrgChains uses r8,r9, so
   mov r5,r9           //   save them once in wrapper function
   push {r4,r5}        //stack contains r8,r9

   ldr curr,[bp,gNxt]  //curr=bp->next
   cmp curr,#0         //if there's only one element,
   bEQ QRET2           //   cleanup and return

   #define LSZ 2*4     //two local variables
   sub SP,SP,#LSZ      //allocate space for UB and ML
   add pML,SP,#0       //save addr of ML
   add pUB,SP,#4       //save addr of UB
   /*
      Stack frame:
      mrgSort(bp,gNxt,gInt)
         [LR,r7-r4,r9-r8]
         [UB]
         [ML]
   */
   #define ASZ 2*4     //two args pushed to stack
   push {pML,pUB}      //store &ML and &UB as arguments for mrgSortR
   /*
      mrgSort(bp,gNxt,gInt)
         [LR,r7-r4,r9-r8]
         [UB]
         [ML]
         [&pUB]
         [&pML]
         mrgSortR(N,bp,gN,gI,&ML,&UB)
   */
   #undef pML          //pML on stack, frees r4
   #undef pUB          //pUB on stack, frees r5
   #define sorted r4   //whether list is sorted or not
   #define prevLen r5  //value of previous node
   movs sorted,#1      //assume sorted to start
   movs N,#1           //init N=1;
   ldr prevLen,[bp,gI] //prevLen = bp->value;
   //we know N>=2, so begin loop
   //loop to count num nodes and check if sorted
   LOOP:
      ldr currLen,[curr,gI]  //currLen = curr->val;
      subs r5,currLen,prevLen//calculate currLen-prevLen, using r5 as mask
      asrs r5,r5,#31         //r5 = (prevLen > currLen) ? 0xFFFFFFFF:0
      mvns r5,r5             //r5 = (prevLen > currLen) ? 0:0xFFFFFFFF
                             //   same as (prevLen <= currLen) ? 0xFFFFFFFF:0
      ands sorted,r5         //sorted &= (prevLen <= currLen) ? 1:0
      movs prevLen,currLen   //prevLen = currLen (update prevLen) 
      adds N,#1              //N++; counts num nodes
      ldr curr,[curr,gNxt]   //curr=curr->next;
      cmp curr,#0            //if(curr)
   bNE LOOP                  //   continue;

   cmp sorted,#0       //if(sorted)
   bNE DONE            //   return bp;
   
   //             r0 r1 r2   r3 SP  SP+4
   //call mrgSortR(N,bp,gNxt,gI,&ML,&UB)
   //current reg config:
   //r0: bp goes to -> r1
   //r1: gNxt goes to -> r2
   //r2: gI goes to -> r3
   //r3: N goes to -> r0
   //swap bp,N: N to r0, bp to r3
   eors r0,r3          //r0 = bp^N
   eors r3,r0          //r3 = N^(bp^N) = (N^N)^bp = 0^bp = bp
   eors r0,r3          //r0 = bp^(bp^N) = (bp^bp)^N = 0^N = N
   //swap gNxt,bp: bp to r1, gNxt to r3
   eors r1,r3          //r1 = gNxt^bp
   eors r3,r1          //r3 = bp^(gNxt^bp) = (bp^bp)^gNxt = 0^gNxt = gNxt
   eors r1,r3          //r1 = gNxt^(gNxt^bp) = (gNxt^gNxt)^bp = 0^bp = bp
   //swap gI,gNxt: gNxt to r2, gI to r3
   eors r2,r3          //r2 = gI^gNxt
   eors r3,r2          //r3 = gNxt^(gI^gNxt) = (gNxt^gNxt)^gI = 0^gI = gI
   eors r2,r3          //r2 = gI^(gI^gNxt) = (gI^gI)^gNxt = 0^gNxt = gNxt
   bl mrgSortR         //mrgSortR(N,bp,gNxt,gI,&ML,&UB)
   
   DONE:
   /*
      mrgSort(bp,gNxt,gInt)
         [LR,r7-r4,r9-r8]
         [UB]
         [ML]
         [&pUB]
         [&pML]
   */
   #define TSZ LSZ+ASZ //space to free = size of local vars + mrgSortR args
   add SP,SP,#TSZ      //deallocate stack vars
   QRET2:
   pop {r4,r5}         //pop caller's r8,r9 from stack
   mov r8,r4           //restore caller's r8
   mov r9,r5           //restore caller's r9
pop {r4-r7,PC}         //restore caller's r4-r7 and return

QRET1:
bx LR

