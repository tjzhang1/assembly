
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

#define FSZ 4*5           //stack holds 5 regs
#define arg4 FSZ          //position on stack of arg4 (ML)
#define arg5 FSZ+4        //position on stack of arg5 (UB)

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
   //r4: pEndL,ML
   //r5: pEndR
   //r6: nr,lb
   //r7: UB
   push {r4-r7,LR}        //store caller's variables
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
   sub SP,SP,#8           //allocate space for endL and endR
   add pEndL,SP,#0        //save address of endL as SP[0]
   add pEndR,SP,#4        //save address of endR as SP[4]
   /*
      Stack frame:
      mrgSortWrapper(...)
      [...]
      [UB]
      [ML]
         mrgSortR(N,bp,gN,gI,ML,UB) (this function)
         [LR,r7-r4,r8]
         [endR]
         [endL]
         [UB]
         [&endL]
            mrgSortR(nl,bp,gN,gI,&endL,UB)
   */
   #define UB r7
   ldr UB,[SP,#arg5]      //load UB
   push {pEndL,UB}        //push args to stack: pEndL=SP[0],UB=SP[4]
   bl mrgSortR            //mrgSortR(nl,bp,gN,gI,&endL,UB)
   #undef nl              //nl no longer needed

   /*
         mrgSortR(N,bp,gN,gI,ML,UB) (this function)
         [LR,r7-r4,r8]
         [endR]
         [endL]
         [UB]
         [&endR]
            mrgSortR(nr,*UB,gN,gI,&endR,UB)
   */
   //need to swap lb (r0) with nr (r6)
   eors r0,nr             //r0 = lb^nr
   eors r6,r0             //r6 = nr^(lb^nr) = (nr^nr)^lb = 0^lb = lb
   #define lb r6          //left buffer
   eors r0,lb             //r0 = lb^(lb^nr) = (lb^lb)^nr = 0^nr = nr
   #undef nr              //r6 now lb, arg0=nr
   str r1,[UB]            //arg1=*UB for mrgSortR
   str pEndR,[SP]         //overwrite SP[0]=&endL with &endR
   #undef UB              //UB no longer needed
   bl mrgSortR            //mrgSortR(nr,*UB,gN,gI,&endR,UB)
   
   /*
         mrgSortR(N,bp,gN,gI,ML,UB) (this function)
         [LR,r7-r4]
         [endR]
         [endL]
         [ML]
            mrgChains(rb,lb,gN,gI,ML,endL,endR)
   */
   movs r1,lb             //arg1=lb for mrgChains
   #undef pEndL           //free r4
   #undef pEndR           //free r5
   add SP,SP,#4           //dealloc &endR
   ldr r4,[SP,#arg4]      //load ML from stack (arg4)
   str r4,[SP]            //overwrite SP[0]=UB with ML
   bl mrgChains           //mrgChains(rb,lb,gN,gI,ML,endL,endR)
   
   //cleanup 
   add SP,SP,#12          //dealloc ML,endL,endR
pop {r4-r7,PC}            //restore caller's r4-r7, and put LR->PC
//bx LR (implicit)

#define ML r4
#define UB r5
#define next r6
BASE:                  //base case: if(N==1)
   ldr UB,[SP,#arg4]      //load UB from stack
   ldr ML,[SP,#arg5]      //load ML from stack
   ldr next,[bp,gNxt]     //load bp->next
   str next,[UB]          //*UB = bp->next; points to next unsorted node
   str bp,[ML]            //*ML = bp; points to last sorted node (bp)
   eors r6,r6             //r6 = #0 = r6^r6
   str next,[bp,gNxt]     //bp->next = NULL; for mrgChains to work
   movs r0,bp             //return bp;
pop {r4-r7,PC}            //restore caller's r4-r7, and put LR->PC
//bx LR (implicit)
