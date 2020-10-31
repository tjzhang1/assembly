
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
   lb = mrgSortRC(leftHalf,bp,gNxt,gI,&endL,&rb);
   //rb set to next unsorted spot
   //endL set to last node in lb
   //mrgSort on right nodes, utilizing the new rb value from prev mrgSort
   rb = mrgSortRC(rightHalf,rb,gNxt,gI,&endR,UB);
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

#define FSZ 0             //stack holds 0 regs
#define arg4 FSZ          //position on stack of arg4 (ML)
#define arg5 FSZ+4        //position on stack of arg5 (UB)

.text
.thumb
.syntax unified
.local BASE
.global mrgSortR

mrgSortR:
   /* Reg configuration A */
   #define N r0           //num nodes to sort
   #define bp r1          //base pointer of input unsorted list
   #define gNxt r2        //num bytes from start of struct for ->next ptr
   #define gI r3          //num bytes from start of struct for ->value

   //bp is always non-NULL and N>0 thanks to wrapper function
   cmp N,#1               //if(N==1)
      bEQ BASE            //   perform base case code
   movs r0,bp             //return head
bx LR

#define UB r4
#define ML r5
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
bx LR
