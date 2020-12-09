
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
   push {r4-r11,LR}    //save caller's variables
   //wrapper only uses r4-r7, but mrgChains uses r8,r9, so
   //   save them once in wrapper function
   ldr curr,[bp,gNxt]  //curr=bp->next
   cmp curr,#0         //if there's only one element,
   bEQ DONE            //   cleanup and return
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
   bNE LOOP                  //   continue;

   cmp sorted,#0       //if(sorted)
   bNE DONE            //   return bp;

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
pop {r4-r11,PC}        //restore caller's r4-r9 and return

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
   else if (N == 2)  //N = 2 case
   {
      word_t *n0=bp, *n1=bp->next; //node 1 and node 2
      int n0Val=n0->len, n1Val=n1->len;

      *UB = n1->next;
      if(n0Val <= n1Val)  //already sorted
      {
         *ML = n1;
         n1->next = NULL;
         return n0;
      }
      //if not sorted, then swap n0 and n1
      n1->next = n0;
      n0->next = NULL;
      *ML = n0;
      return n1;
   }
   //bp is always non-NULL and N>0 thanks to wrapper function
   word_t *lb, *rb;
   word_t *endL, *endR;  //ptrs to last nodes of lb and rb
   int minL, minR;       //holds top (min) val of lb and rb
   int maxL, maxR;       //holds last (max) val of lb and rb
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
   maxL = endL->len;
   minR = rb->len;
   minL = lb->len;
   maxR = endR->len;
   //if max val in lb <= min val in rb, append rb to lb
   if(maxL <= minR)
   {
      *ML = endR;
      endL->next=rb;
      return lb;
   }
   //if max val in rb <= smallest val in lb, append lb to rb
   else if (maxR <= minL)
   {
      *ML = endL;
      endR->next=lb;
      return rb;
   }
   //Make it so lb will always run out of nodes before rb. This always occurs
   //  if max lb <= max rb. Benefits: simplifying mrgChains loops
   //  Since lb runs out of nodes before rb, then only need to check when lb is
   //  empty before appending the rest of rb onto lb.
   //if max lb > max rb, swap lb and rb to make this assumption true
   if (maxL > maxR)
   {
      //last node of merged list is endL now
      *ML = endL;
      //call mrgChains with lb and rb swapped
      return mrgChainsC(rightHalf,leftHalf,rb,lb,minR,minL);
   }
   else
   {
      //the last node of merged list is endR
      *ML = endR;
      //call mrgChains with lb and rb in standard order
      return mrgChainsC(leftHalf,rightHalf,lb,rb,minL,minR);
   }
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
.local BASE1,BASE2,CAT_LB,CAT_RB,ELSE
.global mrgSortR

#define nl r0             //num nodes in left half
#define rb r0             //right buffer
#define bp r1             //base pointer of input unsorted list
#define maxL r1           //value of endL->len
#define gNxt r2           //num bytes from start of struct for ->next ptr
#define gI r3             //num bytes from start of struct for ->value
#define endL r4           //addr of last node of left buffer
#define endR r5           //addr of last node of right buffer
#define nr r6             //num nodes in right half
#define lb r6             //left buffer
#define pUB r7            //ptr to addr of next node after sorted list
#define minL r8           //value of lb->len
#define minR r9           //value of rb->len
#define maxR r10          //value of endR->len
#define N r11             //num nodes to sort
#define pML r12           //ptr to addr of last node in mrged list

mrgSortR:
   /* Reg configuration */
   //r0: nl,nr,lb,rb
   //r1: bp,maxL
   //r2: gNxt
   //r3: gI
   //r4: endL             
   //r5: endR
   //r6: nr,lb
   //r7: pUB              //Don't need to save pUB because it is unchanging
   //r8: minL             //Don't need to save r8-r10 because they don't need
   //r9: minR             //   to maintain value across calls to mrgSortR
   //r10: maxR            
   //r11: N               //num nodes
   //r12: pML             //pML is unchanging so no need to save
   push {r4-r6,r11,LR}    //store caller's variables
   mov N,r0               //store r0 into r11
   //bp is always non-NULL and N>0 thanks to wrapper function
   cmp N,#1               //if(N==1)
   bEQ BASE1              //   perform base case code
   cmp N,#2               //if(N==2)
   bEQ BASE2              //   perform base case for N==2 code
   lsrs nl,N,#1           //nl = N >> 1 = N/2
   subs nr,N,nl           //nr = N-nl = N-N/2
                          //         r0 r1 r2 r3 r12 r7
   bl mrgSortR            //mrgSortR(nl,bp,gN,gI,ML,UB)
   ldr endL,[pML]         //endL = *ML = addr of last node of lb

   //need to swap lb (r0) with nr (r6)
   eors r0,r6             //r0 = lb^nr
   eors r6,r0             //r6 = nr^(lb^nr) = (nr^nr)^lb = 0^lb = lb
   eors r0,r6             //r0 = lb^(lb^nr) = (lb^lb)^nr = 0^nr = nr
   #undef nr              //r6 now lb, arg0=nr
   ldr r1,[pUB]           //arg1=*pUB for mrgSortR
                          //         r0  r1 r2 r3 r12 r7
   bl mrgSortR            //mrgSortR(nr,*UB,gN,gI,ML,UB)
   ldr endR,[pML]         //endR = *ML = addr of last node of rb
   //Current reg configuration:
   //r0: rb
   //r1: ___ (trash) use as maxL
   //r2: gN
   //r3: gInt
   //r4: endL,nL
   //r5: endR,nR
   //r6: lb
   //r7: pUB
   //r8: ___ (trash) use as minL
   //r9: ___ (trash) use as minR
   //r10: __ (trash) use as maxR
   //r11: N
   //r12: pML
   ldr maxL,[endL,gI]     //Val of last node in lb is max
   ldr minR,[rb,gI]       //Val of first node in rb in min
   ldr minL,[lb,gI]       //Val of first node in lb is min
   ldr maxR,[endR,gI]     //Val of last node in rb is max
   cmp maxL,minR          //if (maxL <= minR)
   bLE CAT_RB             //   append rb to lb
   cmp maxR,minL          //if (maxR <= minL)
   bLE CAT_LB             //   append lb to rb
   //Make it so lb will always run out of nodes before rb. This always occurs
   //  if max lb <= max rb. Benefits: simplifies mrgChains loops
   //  Since lb runs out of nodes before rb, then only need to check when lb is
   //  empty before appending the rest of rb onto lb.
   //if max lb > max rb, swap lb and rb to make this assumption true
   #undef nl
   #undef nr
   #define nl r4
   #define nr r5
   cmp maxL,maxR          //if(maxL > maxR)
      bLE ELSE
      str endL,[pML]      //   endL is last node
      lsrs nr,N,#1        //   nr = N/2 (=nl,makes it so r5 holds nl)
      subs nl,N,nr        //   nl = N-N/2 (=nr,makes it so r4 holds nr)
      movs r1,lb          //   set lb to arg1 for mrgChains
      eors minL,minR      //   Swap r8 and r9 so vl in r9 and vr in r8
      eors minR,minL      //   minR = minL^minR^minR = minL (=r9)
      eors minL,minR      //   minL = minR^minL^minL = minR (=r8)
                          //                    r0 r1 r2 r3 r4 r5 r8 r9
      bl mrgChains        //   return mrgChains(rb,lb,gN,gI,nr,nl,vr,vl)
   pop {r4-r6,r11,PC}        //restore caller's r4-r6,r11 and put LR->PC
   ELSE:                  //Otherwise:
      str endR,[pML]      //   endR is last node
      lsrs nl,N,#1        //   nl = N >> 1 = N/2
      subs nr,N,nl        //   nr = N-nl = N-N/2
      movs r1,rb          //   set rb to arg1 for mrgChains
      movs r0,lb          //   set lb to arg0 for mrgChains
                          //                    r0 r1 r2 r3 r4 r5 r8 r9
   bl mrgChains           //   return mrgChains(lb,rb,gN,gI,nl,nr,vl,vr)
pop {r4-r6,r11,PC}        //restore caller's r4-r6,r11 and put LR->PC

CAT_RB:                   //occurs when maxL <= minR, append rb to lb
   str endR,[pML]         //last node of merged list is endR
   str rb,[endL,gNxt]     //append rb to endL
   mov r0,lb              //return lb as first node of merged list
pop {r4-r6,r11,PC}        //restore caller's r4-r6,r11 and put LR->PC

CAT_LB:                   //occurs when maxR <= minL, append lb to rb
   str endL,[pML]         //last node of merged list is endL
   str lb,[endR,gNxt]     //append lb to endR
   //return rb (implicit, already in r0)
pop {r4-r6,r11,PC}        //restore caller's r4-r6,r11 and put LR->PC

#define next r6
BASE1:                    //base case: if(N==1)
   ldr next,[bp,gNxt]     //load bp->next
   str next,[pUB]         //*UB = bp->next; points to next unsorted node
   str bp,[pML]           //*ML = bp; points to last sorted node (bp)
   eors next,next         //next = #0 = next^next
   str next,[bp,gNxt]     //bp->next = NULL; for mrgChains to work
   movs r0,bp             //return bp;
pop {r4-r6,r11,PC}        //restore caller's r4-r6,r11 and put LR->PC

#define n1 r0
BASE2:                    //base case: if(N==2)
   ldr n1,[bp,gNxt]       //load node 1
   ldr minL,[bp,gI]       //load val of node 0
   ldr minR,[n1,gI]       //load val of node 1
   ldr next,[n1,gNxt]     //load n1->next
   str next,[pUB]         //*UB = n1->next; points to next unsorted node
   eors next,next         //next = #0 = next^next
   //If val in bp > val in n1, swap them so they are in order
   cmp minL,minR          //The following IT blocks have the following logic:
      ittee GT            //if(n0 > n1)
      strGT bp,[pML]      //   *ML = n0; last merged node is n0
      strGT next,[bp,gNxt]//   bp->next=NULL;
      strLE n1,[pML]      //   n1->next=bp;
      strLE next,[n1,gNxt]//   return n1;
      ite GT              //else
      strGT bp,[n1,gNxt]  //   *ML = n1;
      movLE r0,bp         //   n1->next=NULL;
                          //   return bp;
pop {r4-r6,r11,PC}        //restore caller's r4-r6,r11 and put LR->PC

#undef rb
#undef bp
#undef gNxt
#undef gI
#undef lb
#undef endL
#undef endR
#undef nl 
#undef nr
#undef pML
#undef pUB
#undef next
#undef maxL
#undef minL
#undef maxR
#undef minR
#undef n1

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
   int nL,                //r4
   int nR,                //r5
   int vl,                //r8
   int vr                 //r9
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
.local CLEANUP,ELSE1,TAKE_L,TAKE_R,CONT_L,CONT_R
.global mrgChains

#define lb r0             //holds the top of left chain
#define rb r1             //holds the top of right chain
#define gNxt r2           //holds offset in bytes to ptr to next node
#define gI r3             //holds offset in bytes to value of node
#define nl r4             //num remaining nodes in lb
#define nr r5             //num remaining nodes in rb
#define head r6           //holds the base ptr of sorted list
#define curr r7           //holds curr node of sorted list
#define vl r8             //val of first node in lb
#define vr r9             //val of first node in rb
#define nextL r10         //preload lb->next
#define nextR r11         //preload rb->next
#define nextLV r12        //preload lb->next->len
#define nextRV r14        //preload rb->next->len

mrgChains:
   //Reg configuration
   //r0: lb
   //r1: rb
   //r2: gNxt
   //r3: gInt
   //r4: nl
   //r5: nr
   //r6: head
   //r7: curr             //needs to be saved, mrgSortR uses r7
   //r8: vl               //doesn't need to be saved, mrgSortR doesn't reuse
   //r9: vr               //doesn't need to be saved, mrgSortR doesn't reuse
   //r10: nextL
   //r11: nextR
   //r12: nextLV          //needs to be saved, mrgSortR uses r12
   //r14: nextRV
   push {r7,r12,LR}       //save caller's variables
   //ASSUMES: lb will run out of nodes before rb; i.e. maxlb <= maxrb
   ldr nextL,[lb,gNxt]    //begin loading lb->next
   ldr nextR,[rb,gNxt]    //begin loading rb->next
   subs nr,#1             //pre-subtract 1 from nr, saving me a cmp nr,#1 later
   //Pre-subtracting from nl unnecessary, since I check nl against both 1
   //   and 0 in my loops
   cmp nl,#1              //if (nl > 1)
      it GT
      ldrGT nextLV,[nextL,gI]//load lb->next->len safely
   //It is impossible for rb to have only one element. If this were true, then
   //   since maxlb <= (maxrb = minrb), mrgSortR would have handled this case.
   //   This means I get to load nextRV without checking if nr is zero.
   ldr nextRV,[nextR,gI]  //load rb->next->len
   //Initialize base pointer head
   cmp vl,vr              //if(init_valLeft <= init_valRight)
   bGT ELSE1
   subs nl,#1             //   nL--
   mov head,lb            //   head = lb; (set head of sorted list)
   mov curr,lb            //   keep track of current pos w/o losing head ptr
   mov lb,nextL           //   update lb to lb->next
   mov vl,nextLV          //   update vl to lb->next->len
   //Don't need to check if lb empty because that would mean lb has only one
   //   element < minR. If this were true, then (minL=maxL) < minR, a case
   //   which was already handled in mrgSortR.
   ldr nextL,[nextL,gNxt] //   preload lb->next for next iteration
   cmp nl,#1              //   if (nl > 1)
      it GT
      ldrGT nextLV,[nextL,gI]//preload lb->next->len safely
   cmp vl,vr              //   if(valLeft <= valRight)
   bLE CONT_L             //      continue appending left nodes
                          //   else append right nodes
   TAKE_R:
      subs nr,#1             //nR--;
      str rb,[curr,gNxt]     //set nxt node in sorted list = rb
      mov curr,rb            //curr = rb; (set head of sorted list)
      mov rb,nextR           //update rb to lb->next
      mov vr,nextRV          //update vr to lb->next->len
      //don't need to check for !rb because lb runs out of nodes first
      ldr nextR,[nextR,gNxt] //preload rb->next for next iteration
         it NE               //if(nR) (if there are nodes left)
         ldrNE nextRV,[nextR,gI]//preload lb->next->len safely
      cmp vl,vr              //if (vl <= vr)
      bLE TAKE_L             //   append left nodes
                             //else continue appending right nodes
      CONT_R:
         subs nr,#1             //nR--;
         //don't need to set curr->next because rb already points to rb->next
         mov curr,rb            //curr = rb;
         mov rb,nextR           //update rb to rb->next
         mov vr,nextRV          //update vr to rb->next->len
         //don't need to check for !rb because lb runs out of nodes first
         ldr nextR,[nextR,gNxt] //preload rb->next for next iteration
            it NE               //if(nR) (if there are nodes left)
            ldrNE nextRV,[nextR,gI]//preload lb->next->len safely
         cmp vl,vr              //if (vl <= vr)
         bLE TAKE_L             //   append left node
      b CONT_R                  //else continue appending right nodes

   ELSE1:                 //else if(init_valLeft > init_valRight)
   subs nr,#1             //   nR--
   mov head,rb            //   head = rb; (set head of sorted list)
   mov curr,rb            //   keep track of current pos w/o losing head ptr
   mov rb,nextR           //   update rb to rb->next
   mov vr,nextRV          //   update vr to rb->next->len
   //don't need to check for !rb because lb runs out of nodes first
   ldr nextR,[nextR,gNxt] //   preload rb->next for next iteration
      it NE               //   if (nR) (if there are nodes left)
      ldrNE nextRV,[nextR,gI]//preload rb->next->len safely
   cmp vl,vr              //   if(valLeft > valRight)
   bGT CONT_R             //      continue appending right nodes
                          //   else append left nodes
   TAKE_L:
      subs nl,#1             //nl--;
      str lb,[curr,gNxt]     //set nxt node in sorted list = lb
      mov curr,lb            //curr = lb;
      bEQ CLEANUP            //if(!lb)   break;
      mov lb,nextL           //update lb to lb->next
      mov vl,nextLV          //update vl to lb->next->len
      ldr nextL,[nextL,gNxt] //preload lb->next for next iteration
      cmp nl,#1              //if (nl > 1)
         it GT
         ldrGT nextLV,[nextL,gI]//preload lb->next->len safely
      cmp vl,vr              //if (vl > vr)
      bGT TAKE_R             //   append right nodes
                             //else continue appending left nodes
      CONT_L:
         subs nl,#1             //nl--;
         //don't need to set curr->next because lb already points to lb->next
         mov curr,lb            //curr = lb;
         bEQ CLEANUP            //if(!lb)   break;
         mov lb,nextL           //update lb to lb->next
         mov vl,nextLV          //update vl to lb->next->len
         ldr nextL,[nextL,gNxt] //preload lb->next for next iteration
         cmp nl,#1              //if (nl > 1)
            it GT
            ldrGT nextLV,[nextL,gI]//preload lb->next->len safely
         cmp vl,vr              //if (vl > vr)
         bGT TAKE_R             //   append right nodes
      b CONT_L                  //else continue appending left nodes

   CLEANUP:               //occurs if(!lb) -> add rb to curr
   str rb,[curr,gNxt]     //append rb to current node
   movs r0,head           //return head
pop {r7,r12,PC}           //restore the rest of caller's variables

