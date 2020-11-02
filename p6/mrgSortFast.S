
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

QRET1:  //occurs when bp=NULL, just return bp
bx LR

#undef bp
#undef gNxt
#undef gI
#undef N
#undef sorted
#undef prevLen
#undef currLen
#undef curr
#undef TSZ
#undef LSZ
#undef ASZ

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
   void *ML,              //SP[0]
   void *UB               //SP[4]
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
#define bp r1             //base pointer of input unsorted list
#define gNxt r2           //num bytes from start of struct for ->next ptr
#define gI r3             //num bytes from start of struct for ->value
#define pEndL r4          //ptr to addr of last node of left buffer
#define pEndR r5          //ptr to addr of last node of right buffer

mrgSortR:
   /* Reg configuration */
   //r0: N,nl,nr,lb,rb
   //r1: bp,*UB
   //r2: gNxt
   //r3: gI
   //r4: pEndL,ML         //OUT: set ML ptr to addr of last node in sorted list
   //r5: pEndR
   //r6: nr,lb
   //r7: UB               //OUT: set ptr to addr of next node after sorted list
   push {r4-r7,LR}        //store caller's variables
   #define FSZ 4*5        //stack holds 5 regs

   //bp is always non-NULL and N>0 thanks to wrapper function
   cmp N,#1               //if(N==1)
      bEQ BASE            //   perform base case code
   //perform splitting
   #define nl r0          //num nodes in left half
   #define nr r6          //num nodes in right half
   lsrs nr,N,#1           //nr = N >> 1 = N/2
   subs nl,N,nr           //nl = N-nr = N-N/2
   #undef N               //N no longer needed
   //word_t *endL, *endR;
   #define LSZ 2*4        //save 2 local variables
   sub SP,SP,#LSZ         //allocate space for endL and endR
   add pEndL,SP,#0        //save address of endL as SP[0]
   add pEndR,SP,#4        //save address of endR as SP[4]

   #define UB r7
   #define arg5 FSZ+LSZ+4 //position on stack of UB; arg5
   ldr UB,[SP,#arg5]      //load UB from stack
   #undef arg5
   push {pEndL,UB}        //push args to stack: pEndL=SP[0],UB=SP[4]
   /*
      Stack frame:
      mrgSortWrapper(...)
      [...]
      [UB]
      [ML]
         mrgSortR(N,bp,gN,gI,ML,UB) (this function)
         [LR,r7-r4] (caller's regs)
         [endR]
         [endL]
         [UB]
         [&endL]
            mrgSortR(nl,bp,gN,gI,&endL,UB)
   */
   bl mrgSortR            //mrgSortR(nl,bp,gN,gI,&endL,UB)
   #undef nl              //nl no longer needed

   //need to swap lb (r0) with nr (r6)
   eors r0,nr             //r0 = lb^nr
   eors r6,r0             //r6 = nr^(lb^nr) = (nr^nr)^lb = 0^lb = lb
   #define lb r6          //left buffer
   eors r0,lb             //r0 = lb^(lb^nr) = (lb^lb)^nr = 0^nr = nr
   #undef nr              //r6 now lb, arg0=nr
   ldr r1,[UB]            //arg1=*UB for mrgSortR
   str pEndR,[SP]         //overwrite SP[0]=&endL with &endR
   #undef UB              //UB no longer needed
   /*
         mrgSortR(N,bp,gN,gI,ML,UB) (this function)
         [LR,r7-r4] (caller's regs)
         [endR]
         [endL]
         [UB]
         [&endR]
            mrgSortR(nr,*UB,gN,gI,&endR,UB)
   */
   bl mrgSortR            //mrgSortR(nr,*UB,gN,gI,&endR,UB)

   movs r1,r0             //set rb to arg1 for mrgChains
   movs r0,lb             //arg0=lb for mrgChains
   #undef pEndL           //free r4
   #undef pEndR           //free r5
   add SP,SP,#8           //dealloc &endR and UB
   #define arg4 FSZ+LSZ   //position on stack of ML; arg4
   ldr r4,[SP,#arg4]      //load ML from stack (arg4)
   #undef arg4
   push {r4}              //push ML to stack
   /*
         mrgSortR(N,bp,gN,gI,ML,UB) (this function)
         [LR,r7-r4]
         [endR]
         [endL]
         [ML]
            mrgChains(lb,rb,gN,gI,ML,endL,endR)
   */
   bl mrgChains           //mrgChains(lb,rb,gN,gI,ML,endL,endR)

   //cleanup
   add SP,SP,#12          //dealloc ML,endL,endR
pop {r4-r7,PC}            //restore caller's r4-r7, and put LR->PC
//bx LR (implicit)

#define ML r4
#define UB r5
#define next r6
#define arg4 FSZ          //position on stack of ML; arg4
#define arg5 FSZ+4        //position on stack of UB; arg5
BASE:                  //base case: if(N==1)
   ldr ML,[SP,#arg4]      //load ML from stack
   ldr UB,[SP,#arg5]      //load UB from stack
   ldr next,[bp,gNxt]     //load bp->next
   str next,[UB]          //*UB = bp->next; points to next unsorted node
   str bp,[ML]            //*ML = bp; points to last sorted node (bp)
   eors next,next         //next = #0 = next^next
   str next,[bp,gNxt]     //bp->next = NULL; for mrgChains to work
   movs r0,bp             //return bp;
pop {r4-r7,PC}            //restore caller's r4-r7, and put LR->PC
//bx LR (implicit)

#undef bp
#undef gNxt
#undef gI
#undef lb
#undef ML
#undef UB
#undef next
#undef arg4
#undef arg5
#undef FSZ
#undef LSZ

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
   void *ML,              //SP[0]
   void *endL,            //SP[4]
   void *endR             //SP[8]
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
#define vl r4          //holds first val of lb = min val of lb
#define maxR r5        //holds val of endR = max val of rb
#define endL r6        //holds last node of lb
#define endR r7        //holds last node of rb

#define FSZ 4*4        //4 regs stored
#define arg4 FSZ       //offset of first arg on SP (ML)
#define arg5 FSZ+4     //offset of 2nd arg on SP (endL)
#define arg6 FSZ+8     //offset of 3rd arg on SP (endR)

mrgChains:
   //Reg configuration A
   //r0: lb
   //r1: rb
   //r2: gNxt
   //r3: gInt
   //r4: vl / vr
   //r5: maxR / maxL
   //r6: endL
   //r7: endR
   push {r4-r7}           //save caller's variables
   ldr endL,[SP,#arg5]    //load endL argument from stack
   ldr endR,[SP,#arg6]    //load endR argument from stack
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
      bLE QRET3           //   concatenate lb to rb and return
   mov r8,vl              //save vl by moving to higher reg
   #undef vl              //reassign vl for vr
   #undef maxR            //free r5 for maxL
   #define vr r4          //holds first val of rb = min val of rb
   #define maxL r5        //holds val of endL = max val of lb
   ldr vr,[rb,gInt]       //load int value from rb
   ldr maxL,[endL,gInt]   //load int value from endL
   cmp maxL,vr            //if(max(lb) <= min(rb))
      bLE QRET4           //   concatenate rb to lb and return
   mov r5,r8              //restore vl from high reg
   mov r8,endL            //store endL to high reg
   mov r9,endR            //store endR to high reg
   movs r7,lb             //make room for new list head
   #undef maxL            //free r5 for vl
   #undef endL            //free r6 for curr
   #undef endR            //free r7 for lb
   #undef lb              //reassign lb
   #define head r0        //holds the base ptr of sorted list
   #define vl r5          //new spot for vl
   #define curr r6        //holds curr node of sorted list
   #define lb r7          //new spot for lb
   #define endL r8        //new spot for endL
   #define endR r9        //new spot for endR

   //Reg configuration B
   //r0: head
   //r1: rb
   //r2: gNxt
   //r3: gInt
   //r4: vr
   //r5: vl
   //r6: curr
   //r7: lb
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
#undef vr
#undef vl
#undef curr
#undef lb
#undef endL
#undef endR

/* Reg configuration A */
#define lb r0             //holds the top of left chain
#define rb r1             //holds the top of right chain
#define gNxt r2           //holds offset in bytes to ptr to next node
#define ML r4             //updated to point to last node
#define endL r6           //holds last node of lb
#define endR r7           //holds last node of rb
QRET3:                 //occurs if(max(rb) <= min(lb)) -> add lb to rb
   ldr ML,[SP,#arg4]      //load ML from stack
   str endL,[ML]          //*ML = endL
   str lb,[endR,gNxt]     //append lb to end of rb
   movs r0,rb             //return rb
   pop {r4-r7}            //restore the rest of caller's variables
bx LR

QRET4:                 //occurs if(max(lb) <= min(rb)) -> add rb to lb
   ldr ML,[SP,#arg4]      //load ML from stack
   str endR,[ML]          //*ML = endL
   str rb,[endL,gNxt]     //append rb to end of lb
   pop {r4-r7}            //restore the rest of caller's variables
bx LR                     //return lb
#undef lb
#undef endL
#undef endR

/* Reg configuration B */
#define curr r6           //holds curr node of sorted list
#define lb r7             //holds the top of left chain
#define endL r8           //holds last node of lb
#define endR r9           //holds last node of rb
DONE1:                 //occurs if(!lb) -> add rb to curr
   ldr ML,[SP,#arg4]      //load ML from stack
   mov r5,endR            //mov endR from higher register
   str r5,[ML]            //*ML = endR
   str rb,[curr,gNxt]     //append rb to current node
   pop {r4-r7}            //restore the rest of caller's variables
bx LR                     //return head (in r0 already)

DONE2:                 //occurs if(!rb) -> add lb to curr
   ldr ML,[SP,#arg4]      //load ML from stack
   mov r5,endL            //mov endL from higher register
   str r5,[ML]            //*ML = endL
   str lb,[curr,gNxt]     //append lb to current node
   pop {r4-r7}            //restore the rest of caller's variables
bx LR                     //return head (in r0 already)
