#define bp r0                //holds ptr to head as well as current node
#define pMinPrev r1          //holds ptr to ptr to node before minimum
#define nextgap r2           //holds num bytes from start of node to next ptr
#define igap r3              //holds num bytes from start of node to integer
#define minPrev r12          //holds ptr to node before minimum
#define minNode r4           //holds minimum node
#define minVal r5            //holds minimum value
#define val r6               //holds val of current node
#define prev r7              //holds ptr to node before current node

/*
r0                         r0            r1                    r2
void *findIntMinNode(void *bp,void *minPrev,unsigned long nextgap,
                   r3
   unsigned long igap);
{
   struct someStruct *minNode, *prev, *prevMin;
   int minval;
   if(minPrev)
      *minPrev = NULL;
   if(!bp)
      return NULL;
   
   minNode = bp;
   //minVal = bp->int_val;
   //in ARM will look something like:
   //minVal = *((void *)((unsigned long)bp + igap));
   minVal = bp->int_val;
   prev = bp;
   //bp = bp->next;
   //in ARM will look something like:
   //bp = *((void *)((unsigned long)bp + nextgap));
   for(bp=bp->next; bp; bp=bp->next))
   {
      int val = bp->int_val;
      if(val < minVal)
      {
         minVal = val;
         minNode=bp;
         prevMin=prev;
      }
      prev=bp;
   }
   if(minPrev)
      *minPrev=prevMin;
   return minNode;
}

This function finds the node in the linked list starting with bp with minimum
int value. The parameter nextgap is the offset in bytes from the start of the
struct that holds the pointer to the next node. Similarly, igap is the offset
in bytes that is used to find the integer value.
RETURNS: the address of the minimum node. If minPrev isn't NULL, we point this
 to the address of the node before the minimum node.
*/
.text 
.arm
.local LOOP, END
.globl findIntMinNode
findIntMinNode:
   //prologue
   eor minPrev,minPrev,minPrev  //init minPrev to NULL i.e. set to 0
   cmp pMinPrev,#0              //if(pMinPrev)
   strNE minPrev,[pMinPrev]     //   *pMinPrev=NULL;
   cmp bp,#0                    //if(!bp)
   bxEQ lr                      //   return NULL; (since bp=r0=NULL)
   stmDB SP!,{r4-r7}            //save caller's r4-r7

   //body
   mov minNode,bp               //init minNode to first node
   ldr minVal,[bp,igap]         //init minVal to the int in first node
   mov prev,bp                  //init prev=bp
   ldr bp,[bp,nextgap]          //bp=bp->next
   cmp bp,#0                    //if(!bp)
   bEQ END                      //   skip loop
/*
   for(bp=bp->next; bp; bp=bp->next))
   {
      int val = bp->int_val;
      if(val < minVal)
      {
         minVal = val;
         minNode=bp;
         prevMin=prev;
      }
      prev=bp;
   }
Loop through each node and compare its value to the minVal. Updave minVal,
minNode, and prevMin accordingly.
*/
   LOOP:
      ldr val,[bp,igap]         //load the value of current node
      cmp val,minVal            //if current value < minVal
      movLT minVal,val          //   minVal = current value
      movLT minNode,bp          //   minNode = current node
      movLT minPrev,prev        //   minPrev = prev
      mov prev,bp               //update prev to current node
      ldr bp,[bp,nextgap]       //bp=bp->next (increment current node)
      cmp bp,#0                 //if(bp!=0)
   bNE LOOP                     //   continue loop
   cmp pMinPrev,#0              //if(pMinPrev)
   strNE minPrev,[pMinPrev]     //   *pMinPrev=minPrev;
   END:  //code execution jumps here if only one node
   mov r0,minNode               //return minNode

   //epilogue
   ldmIA SP!,{r4-r7}            //restore caller's r4-r7
bx LR
