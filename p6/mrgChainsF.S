#define head r0           //holds the head ptr of sorted list
#define rb r1             //holds the top of right chain
#define gNxt r2           //holds offset in bytes to ptr to next node
#define gInt r3           //holds offset in bytes to value of node
#define lb r12            //holds the top of left chain
#define vl r4             //holds val of lb (vl=valLeft)
#define vr r5             //holds val of rb (vr=valRight)
#define curr r6           //holds curr node of sorted list

/*
r0                    r0        r1       r2         r3
void *mrgChains
(
   void *lb,
   void *rb,
   uint gNxt,
   uint gInt,
   
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
.arm
.local DONE,LOOP,ELSE1,ELSE2
.global mrgChains

mrgChains:
   //prologue:
   stmDB SP!,{r4-r6}      //save caller's variables
   mov lb,r0              //save r0 to make room for head
   
   //body:
   ldr vl,[lb,gInt]       //load int value from lb
   ldr vr,[rb,gInt]       //load int value from rb
   cmp vl,vr              //if(valLeft <= valRight)
      bGT ELSE1
      mov head,lb         //   head = lb; (set head of sorted list)
      ldr lb,[lb,gNxt]    //   lb=lb->next; (increment lb)
      mov curr,head       //   keep track of current pos w/o losing head ptr
      cmp lb,#0           //   if(!lb)
      bEQ DONE            //      skip loop
      ldr vl,[lb,gInt]    //   else update int value from lb
      b LOOP              //   begin loop
   ELSE1:                 //else if(valLeft > valRight)
      mov head,rb         //   head = rb; (set head of sorted list)
      ldr rb,[rb,gNxt]    //   rb=rb->next; (increment rb)
      mov curr,head       //   keep track of current pos w/o losing head ptr
      cmp rb,#0           //   if(!rb)
      bEQ DONE            //      skip loop
      ldr vr,[rb,gInt]    //   else update int value from rb, begin loop

   LOOP:
      cmp vl,vr              //if(valLeft <= valRight)
         bGT ELSE2
         str lb,[curr,gNxt]  //   set nxt node in sorted list = lb 
         mov curr,lb         //   curr = lb; (curr=curr->next)
         ldr lb,[lb,gNxt]    //   lb=lb->next;
         cmp lb,#0           //   if(!lb)
         bEQ DONE            //      break;
         ldr vl,[lb,gInt]    //   else update int value from lb
         b LOOP              //      continue;
      ELSE2:                 //else if(valLeft > valRight)
         str rb,[curr,gNxt]  //   set nxt node in sorted list = rb
         mov curr,rb         //   curr = rb; (curr=curr->next)
         ldr rb,[rb,gNxt]    //   rb=rb->next;
         cmp rb,#0           //   if(!rb)
         bEQ DONE            //      break;
         ldr vr,[rb,gInt]    //   else update int value from rb
   b LOOP                    //      continue;

   DONE:
   //cleanup; at this pt lb=0 or rb=0, but not both
   //append remaining nodes to the end of curr
   //since they are mutually exclusive, bitwise OR them together
   orr lb,rb              //lb = (lb) ? lb:rb
   str lb,[curr,gNxt]     //curr->next = (lb) ? lb:rb
   //epilogue:
   ldmIA SP!,{r4-r6}      //restore caller's variables
bx LR
