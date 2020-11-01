
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
   //r4: pEndL,ML
   //r5: pEndR
   //r6: nr,lb
   //r7: UB
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
