
//Table of contents:
// mrgSort: around line 10
// mrgSortR: around line 170
// mrgChains: around line 290 

//-----------------------------------------------------------------------------
//                   mrgSort(bp,gN,gI)
//-----------------------------------------------------------------------------

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
.local LOOP,DONE
.global mrgSort

#define bp r0
#define gNxt r1
#define gI r2
#define N r3
#define sorted r4   //whether list is sorted or not
#define prevLen r5  //value of previous node
#define currLen r6
#define curr r7
#define pUB r7
#define pML r12

mrgSort:
   cmp bp,#0           //if bp NULL
      it EQ
      bxEQ LR          //   return bp
   push {r4-r9,LR}     //save caller's variables
   //wrapper only uses r4-r7, but mrgChains uses r8,r9, so
   //   save them once in wrapper function
   ldr curr,[bp,gNxt]  //curr=bp->next
   cmp curr,#0         //if there's only one element,
      it EQ
      bEQ DONE         //   cleanup and return
   mov sorted,#1       //assume sorted to start
   mov N,#1            //start N=1 since there's at least 1 node
   ldr prevLen,[bp,gI] //prevLen = bp->value;
   //we know N>=2, so begin loop
   //loop to count num nodes and check if sorted
   LOOP:
      ldr currLen,[curr,gI]  //currLen = curr->val;
      ldr curr,[curr,gNxt]   //curr=curr->next; begin ld before end of loop
      subs r5,currLen,prevLen//calculate currLen-prevLen, using r5 as mask
      asrs r5,r5,#31         //r5 = (prevLen > currLen) ? 0xFFFFFFFF:0
      mvns r5,r5             //r5 = (prevLen > currLen) ? 0:0xFFFFFFFF
                             //   same as (prevLen <= currLen) ? 0xFFFFFFFF:0
      ands sorted,r5         //sorted &= (prevLen <= currLen) ? 1:0
      mov prevLen,currLen    //prevLen = currLen (update prevLen)
      adds N,#1              //N++; counts num nodes
      cmp curr,#0            //if(curr)
         it NE
   bNE LOOP                  //   continue;

   cmp sorted,#0       //if(sorted)
      it NE
      bNE DONE         //   return bp;

   //             r0 r1 r2   r3 r12 r7
   //call mrgSortR(N,bp,gNxt,gI,&ML,&UB)
   #define LSZ 2*4     //two local variables
   sub SP,SP,#LSZ      //allocate space for UB and ML
   add pML,SP,#0       //save addr of ML
   add pUB,SP,#4       //save addr of UB to r7
   /*
      Stack frame:
      mrgSort(bp,gNxt,gInt)
         [LR,r7-r4,r9-r8]
         [UB]
         [ML]
         mrgSortR(N,bp,gN,gI,&ML,&UB)
   */
   //moving regs around:
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
   //                  //        r0 r1 r2   r3 r12 r7
   bl mrgSortR         //mrgSortR(N,bp,gNxt,gI,&ML,&UB)

   /*
      mrgSort(bp,gNxt,gInt)
         [LR,r7-r4,r9-r8]
         [UB]
         [ML]
   */
   add SP,SP,#LSZ      //deallocate local vars
   DONE:
pop {r4-r9,PC}         //restore caller's r4-r9 and return

#undef bp
#undef gNxt
#undef gI
#undef N
#undef pML          //pML on stack, frees r4
#undef pUB          //pUB on stack, frees r5
#undef sorted
#undef prevLen
#undef currLen
#undef curr
#undef LSZ

//-----------------------------------------------------------------------------
//                   mrgSortR(N,bp,gN,gI,ML,UB)
//-----------------------------------------------------------------------------

/*
r0
void *mrgSortRC
(
   unsigned int N,        //r0
   void *bp,              //r1
   unsigned long gNxt,    //r2
   unsigned long gI,      //r3
   void *ML,              //r12
   void *UB               //r7
)
{
   if(N == 1) //base case, sorted
   {
      //increment *UB so it is always pointing to the next unsorted pos
      *UB = bp->next;
      //last node in list is only node in list, bp
      *ML = bp;
      //null terminate
      bp->next = NULL;
      return bp;
   }
   //bp is always non-NULL and N>0 thanks to wrapper function
   word_t *lb, *rb;
   word_t *endL, *endR;  //ptrs to last nodes of lb and rb
   //initialize lengths of left and right buffers for splitting
   unsigned int leftHalf=(N>>1), rightHalf = N-leftHalf;
   //mrgSort on left nodes
   lb = mrgSortRC(leftHalf,bp,gNxt,gI,&endL,UB);
   //rb set to next unsorted spot
   //endL set to last node in lb
   //mrgSort on right nodes, utilizing the new rb value from prev mrgSort
   rb = mrgSortRC(rightHalf,*UB,gNxt,gI,&endR,UB);
   //UB updated to next unsorted spot
   //endR set to last node in rb
   //return mergeLists on left and right
   return mrgChains(lb,rb,gNxt,gI,endL,endR,ML);
}
*/

/*
 * Recursive function performing mrgSort on a linked list of length N
 * ASSUMES: ML and UB non-NULL
 * RETURNS: base pointer to sorted list
 */


.text
.thumb
.syntax unified
.local BASE
.global mrgSortR

#define N r0              //num nodes to sort
#define nl r0             //num nodes in left half
#define bp r1             //base pointer of input unsorted list
#define gNxt r2           //num bytes from start of struct for ->next ptr
#define gI r3             //num bytes from start of struct for ->value
#define endL r4           //addr of last node of left buffer
#define endR r5           //addr of last node of right buffer
#define pML r5            //ptr to addr of last node in mrged list, same as r12
#define nr r6             //num nodes in right half
#define lb r6             //left buffer
#define pUB r7            //ptr to addr of next node after sorted list

mrgSortR:
   /* Reg configuration */
   //r0: N,nl,nr,lb,rb
   //r1: bp,*UB
   //r2: gNxt
   //r3: gI
   //r4: endL             
   //r5: endR,pML         //pML temporarily used here to perform ldr/str
   //r6: nr,lb
   //r7: pUB              //don't need to save pUB because it is unchanging
   //r12: pML             //similarly, pML is unchanging so no need to save
   push {r4-r6,LR}        //store caller's variables

   //bp is always non-NULL and N>0 thanks to wrapper function
   mov pML,r12            //load pML from high reg to do load
   cmp N,#1               //if(N==1)
   bEQ BASE               //   perform base case code
   lsrs nr,N,#1           //nr = N >> 1 = N/2
   subs nl,N,nr           //nl = N-nr = N-N/2
   #undef N               //N no longer needed
                          //         r0 r1 r2 r3 r12 r7
   bl mrgSortR            //mrgSortR(nl,bp,gN,gI,ML,UB)
   ldr endL,[pML]         //endL = *ML = addr of last node of lb
   #undef nl              //nl no longer needed

   //need to swap lb (r0) with nr (r6)
   eors r0,r6             //r0 = lb^nr
   eors r6,r0             //r6 = nr^(lb^nr) = (nr^nr)^lb = 0^lb = lb
   eors r0,r6             //r0 = lb^(lb^nr) = (lb^lb)^nr = 0^nr = nr
   #undef nr              //r6 now lb, arg0=nr
   ldr r1,[pUB]           //arg1=*pUB for mrgSortR
                          //         r0  r1 r2 r3 r12 r7
   bl mrgSortR            //mrgSortR(nr,*UB,gN,gI,ML,UB)
   ldr endR,[pML]         //endR = *ML = addr of last node of rb

   movs r1,r0             //set rb to arg1 for mrgChains
   movs r0,lb             //arg0=lb for mrgChains
                          //                 r0 r1 r2 r3 r12 r4  r5
   bl mrgChains           //return mrgChains(lb,rb,gN,gI,ML,endL,endR)
pop {r4-r6,PC}            //restore caller's r4-r6, and put LR->PC
//bx LR (implicit)

#define next r6
BASE:                     //base case: if(N==1)
   ldr next,[bp,gNxt]     //load bp->next
   str next,[pUB]         //*UB = bp->next; points to next unsorted node
   str bp,[pML]           //*ML = bp; points to last sorted node (bp)
   eors next,next         //next = #0 = next^next
   str next,[bp,gNxt]     //bp->next = NULL; for mrgChains to work
   movs r0,bp             //return bp;
pop {r4-r6,PC}            //restore caller's r4-r7, and put LR->PC
//bx LR (implicit)

#undef bp
#undef gNxt
#undef gI
#undef lb
#undef endL
#undef endR
#undef pML
#undef pUB
#undef next

//-----------------------------------------------------------------------------
//                   mrgChains(lb,rb,gN,gI,ML,endL,endR)
//-----------------------------------------------------------------------------
/*
r0
void *mrgChains
(
   void *lb,              //r0
   void *rb,              //r1
   uint gNxt,             //r2
   uint gInt,             //r3
   void *ML,              //r12
   void *endL,            //r4
   void *endR             //r5
)
{
   void *head, *curr;
   int valLeft, valRight;

   if(valLeft < valRight)
   {
      head = lb;
      lb = lb->next;
   }
   else
   {
      head = rb;
      rb = rb->next;
   }
   curr = head;
   while(lb && rb)
   {
      if(valLeft < valRight)
      {
         curr->next = lb;
         lb = lb->next;
      }
      else
      {
         curr->next = rb;
         rb = rb->next;
      }
      curr = curr->next;
   }
   if(lb)
      curr->next=lb;
   else
      curr->next=rb;
   return head;
}
Performs a sorted-list merge on chains lb and rb.
RETURNS: base of merged sort-order list, assuming lb and rb sorted
ASSUMES: input lists have at least one node, and NULL terminated
*/

.text
.thumb
.syntax unified
.local DONE1,DONE2,LOOP2,ELSE1,ELSE2,QRET3,QRET4
.global mrgChains

#define lb r0          //holds the top of left chain
#define rb r1          //holds the top of right chain
#define gNxt r2        //holds offset in bytes to ptr to next node
#define gInt r3        //holds offset in bytes to value of node
#define endL r4        //holds last node of lb
#define endR r5        //holds last node of rb
#define vl r6          //holds first val of lb = min val of lb
#define vr r6          //holds first val of rb = min val of rb
#define maxL r7        //holds val of endL = max val of lb
#define maxR r7        //holds val of endR = max val of rb
#define FSZ 4*4        //4 regs stored

mrgChains:
   //Reg configuration A
   //r0: lb
   //r1: rb
   //r2: gNxt
   //r3: gInt
   //r4: endL          //doesn't need to be saved, mrgSortR doesn't reuse
   //r5: endR          //doesn't need to be saved, mrgSortR doesn't reuse
   //r6: vl / vr
   //r7: maxR / maxL
   //r12: pML
   push {r6,r7}           //save caller's variables
   /*
      if(maxL <= rbVal)
      {
         *ML = endR;
         endL->next=rb;
         return lb;
      }
      else if (maxR <= lbVal)
      {
         *ML = endL;
         endR->next=lb;
         return rb;
      }
   */
   ldr vl,[lb,gInt]       //load int value from lb
   ldr maxR,[endR,gInt]   //load int value from endR
   cmp maxR,vl            //if(max(rb) <= min(lb))
   bLE QRET3              //   concatenate lb to rb and return
   mov r8,vl              //save vl by moving to higher reg
   #undef vl              //vl now in r8
   #undef maxR            //free r7 for maxL
   ldr vr,[rb,gInt]       //load int value from rb
   ldr maxL,[endL,gInt]   //load int value from endL
   cmp maxL,vr            //if(max(lb) <= min(rb))
   bLE QRET4              //   concatenate rb to lb and return
   //mov regs into better configuration
   mov r7,r8              //restore vl from high reg
   mov r8,endL            //store endL to high reg
   mov r9,endR            //store endR to high reg
   movs r4,lb             //make room for new list head
   #undef maxL            //maxL no longer needed
   #undef endL            //endL now in r8
   #undef endR            //endR now in r9
   #undef lb              //reassign lb to r4
   #define head r0        //holds the base ptr of sorted list
   #define lb r4          //new spot for lb
   #define curr r5        //holds curr node of sorted list
   #define endL r8        //new spot for endL
   #define vl r7          //new spot for vl
   #define endR r9        //new spot for endR

   //Reg configuration B
   //r0: head
   //r1: rb
   //r2: gNxt
   //r3: gInt
   //r4: lb
   //r5: curr
   //r6: vr
   //r7: vl
   //r8: endL
   //r9: endR
   /*
      Initialize base pointer head
   */
   cmp vl,vr              //if(valLeft <= valRight)
      bGT ELSE1
      movs head,lb        //   head = lb; (set head of sorted list)
      ldr lb,[lb,gNxt]    //   lb=lb->next; (increment lb)
      movs curr,head      //   keep track of current pos w/o losing head ptr
      cmp lb,#0           //   if(!lb)
      bEQ DONE1           //      skip loop
      ldr vl,[lb,gInt]    //   else update int value from lb
      b LOOP2             //   begin loop
   ELSE1:                 //else if(valLeft > valRight)
      movs head,rb        //   head = rb; (set head of sorted list)
      ldr rb,[rb,gNxt]    //   rb=rb->next; (increment rb)
      movs curr,head      //   keep track of current pos w/o losing head ptr
      cmp rb,#0           //   if(!rb)
      bEQ DONE2           //      skip loop
      ldr vr,[rb,gInt]    //   else update int value from rb, begin loop

   LOOP2:
      cmp vl,vr              //if(valLeft <= valRight)
         bGT ELSE2
         str lb,[curr,gNxt]  //   set nxt node in sorted list = lb
         movs curr,lb        //   curr = lb; (curr=curr->next)
         ldr lb,[lb,gNxt]    //   lb=lb->next;
         cmp lb,#0           //   if(!lb)
         bEQ DONE1           //      break;
         ldr vl,[lb,gInt]    //   else update int value from lb
         b LOOP2             //      continue;
      ELSE2:                 //else if(valLeft > valRight)
         str rb,[curr,gNxt]  //   set nxt node in sorted list = rb
         movs curr,rb        //   curr = rb; (curr=curr->next)
         ldr rb,[rb,gNxt]    //   rb=rb->next;
         cmp rb,#0           //   if(!rb)
         bEQ DONE2           //      break;
         ldr vr,[rb,gInt]    //   else update int value from rb
b LOOP2                      //      continue;


//Clear reg configuration
#undef head
#undef rb
#undef gNxt
#undef gInt
#undef curr
#undef lb
#undef vr
#undef vl
#undef endL
#undef endR


#define pML r6            //updated to point to last node
/* Reg configuration A */
#define lb r0             //holds the top of left chain
#define rb r1             //holds the top of right chain
#define gNxt r2           //holds offset in bytes to ptr to next node
#define endL r4           //holds last node of lb
#define endR r5           //holds last node of rb
QRET3:                    //occurs if(max(rb) <= min(lb)) -> add lb to rb
   mov pML,r12            //move pML to lower reg for str
   str endL,[pML]         //*ML = endL
   str lb,[endR,gNxt]     //append lb to end of rb
   movs r0,rb             //return rb
   pop {r6,r7}            //restore the rest of caller's variables
bx LR

QRET4:                    //occurs if(max(lb) <= min(rb)) -> add rb to lb
   mov pML,r12            //move pML to lower reg for str
   str endR,[pML]         //*ML = endR
   str rb,[endL,gNxt]     //append rb to end of lb
   pop {r6,r7}            //restore the rest of caller's variables
bx LR                     //return lb
#undef lb
#undef endL
#undef endR

/* Reg configuration B */
#define lb r4             //holds the top of left chain
#define curr r5           //holds curr node of sorted list
#define endL r8           //holds last node of lb
#define endR r9           //holds last node of rb
DONE1:                    //occurs if(!lb) -> add rb to curr
   mov pML,r12            //move pML to lower reg for str
   mov r7,endR            //mov endR from higher register
   str r7,[pML]           //*ML = endR
   str rb,[curr,gNxt]     //append rb to current node
   pop {r6,r7}            //restore the rest of caller's variables
bx LR                     //return head (in r0 already)

DONE2:                    //occurs if(!rb) -> add lb to curr
   mov pML,r12            //move pML to lower reg for str
   mov r7,endL            //mov endL from higher register
   str r7,[pML]           //*ML = endL
   str lb,[curr,gNxt]     //append lb to current node
   pop {r6,r7}            //restore the rest of caller's variables
bx LR                     //return head (in r0 already)
