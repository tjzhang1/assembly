#define head r0           //holds the head ptr of sorted list
#define rb r1             //holds the top of right chain
#define gNxt r2           //holds offset in bytes to ptr to next node
#define gInt r3           //holds offset in bytes to value of node
#define vl r4             //holds val of lb (vl=valLeft)
#define vr r5             //holds val of rb (vr=valRight)
#define curr r6           //holds curr node of sorted list
#define lb r7             //holds the top of left chain

/*
r0                    r0        r1       r2         r3
void *mrgChains(void *lb, void *rb, uint gNxt, uint gInt)
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
.local DONE,LOOP
.global mrgChains

mrgChains:
   //prologue:
   push {r4-r7}           //save caller's variables
   mov lb,r0              //save r0 to make room for head
   
   //body:
   ldr vl,[lb,gInt]       //load int value from lb
   ldr vr,[rb,gInt]       //load int value from rb
   //initialize head ptr
   cmp vl,vr              
      ittee LE            //if(valLeft <= valRight)
      movLE head,lb       //   head = lb; (set head of sorted list)
      ldrLE lb,[lb,gNxt]  //   lb=lb->next;
      movGT head,rb       //else head = rb; (set head of sorted list)
      ldrGT rb,[rb,gNxt]  //   rb=rb->next;
   mov curr,head          //keep track of current pos w/o losing head ptr
   //check to make sure both lists have nodes left
   cmp lb,#0              
      ite NE              //if(!lb) skip loop;
      cmpNE rb,#0         //else if(!rb)
      bEQ DONE            //   skip loop
   //loop to compare nodes from lb and rb and append to sorted ptr
   LOOP:
      ldr vl,[lb,gInt]       //update int value from lb
      ldr vr,[rb,gInt]       //update int value from rb
      cmp vl,vr              
         ittee LE            //if(valLeft <= valRight)
         strLE lb,[curr,gNxt]//   curr->next = lb;
         ldrLE lb,[lb,gNxt]  //   lb=lb->next;
         strGT rb,[curr,gNxt]//else curr->next = rb;
         ldrGT rb,[rb,gNxt]  //   rb=rb->next;
      ldr curr,[curr,gNxt]   //curr=curr->next;
      cmp lb,#0              
         itt NE              //if(lb)
         cmpNE rb,#0         //   && rb)
   bNE LOOP                  //   continue loop

   DONE:
   //cleanup; at this pt lb=0 or rb=0, but not both
   //append remaining nodes to the end of curr
   //since they are mutually exclusive, bitwise OR them together
   orrs lb,rb             //lb = (lb) ? lb:rb
   str lb,[curr,gNxt]     //curr->next = (lb) ? lb:rb
   //epilogue:
   pop {r4-r7}            //restore caller's variables
bx LR
