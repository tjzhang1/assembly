#define sum r0               //holds sum
#define gapNxt r1            //holds offset to next ptr
#define bp r2                //holds ptr to linked list

/*
r0                 r0                r1
int cntNodes(void *bp, unsigned long gapNxt)
{
   int sum=0;
   if(!bp)
      return 0;
   while(bp)
   {
      sum++;
      bp=(void *)((unsigned long)bp + gapNxt);  //i.e. bp=bp->next;
   }
   return sum;
}

This function counts the number of nodes in the linked list starting with bp.
The parameter gapNxt is the offset in bytes from the start of the struct
that holds the pointer to the next node.
RETURNS: number of nodes in the linked linked list. If bp is NULL, returns 0
*/

.text          
.arm            
.local LOOP
.globl cntNodes  
cntNodes:
   mov bp,r0                 //make space for sum
   eor sum,sum,sum           //init sum=0
   cmp bp,#0                 //if(!bp)
   bxEQ lr                   //   return 0

   LOOP:
      add sum,sum,#1         //sum++
      ldr bp,[bp,gapNxt]     //bp=bp->next
      cmp bp,#0              //if(bp!=0)
   bNE LOOP                  //   continue loop
bx lr                        //return sum
