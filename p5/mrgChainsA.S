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
void *mrgChains(void *lb, void *rb, uint gNxt, uint gInt)
{
   //compare first element in order to set newBase
   //set newbase, increment lb or rb
   //loop {
   //check if lb exists
   //   if lb doesn't exist, append rb
   //   break
   //check if rb exists
   //   if rb doesn't exist, append lb
   //   break
   //(lb and rb can't both =NULL because we always append a node one at a time)
   //if both lb and rb exist
   //   compare their values
   //   append lb or rb
   //   increment lb or rb
   //}
   //return base pointer
}
Performs a sorted-list merge on chains lb and rb.
RETURNS: base of merged sort-order list, assuming lb and rb sorted
ASSUMES: input lists have at least one node, and NULL terminated
*/

.text
.arm
.local DONE,LOOP
.global mrgChains

mrgChains:
   //prologue:
   stmDB SP!,{r4-r6}      //save caller's variables
   mov lb,r0              //save r0 to make room for head
   
   //body:
   ldr vl,[lb,gInt]       //load int value from lb
   ldr vr,[rb,gInt]       //load int value from rb
   cmp vl,vr              //if(valLeft <= valRight)
   movLE head,lb          //   head = lb; (set head of sorted list)
   ldrLE lb,[lb,gNxt]     //   lb=lb->next;
   movGT head,rb          //else head = rb; (set head of sorted list)
   ldrGT rb,[rb,gNxt]     //   rb=rb->next;

   mov curr,head          //keep track of current pos w/o losing head ptr
   cmp lb,#0
   cmpNE rb,#0
   bEQ DONE               //   skip loop
   LOOP:
      ldr vl,[lb,gInt]       //update int value from lb
      ldr vr,[rb,gInt]       //update int value from rb
      cmp vl,vr              //if(valLeft <= valRight)
      strLE lb,[curr,gNxt]   //   curr->next = lb;
      ldrLE lb,[lb,gNxt]     //   lb=lb->next;
      strGT rb,[curr,gNxt]   //else curr->next = rb;
      ldrGT rb,[rb,gNxt]     //   rb=rb->next;
      ldr curr,[curr,gNxt]   //curr=curr->next;
      cmp lb,#0
      cmpNE rb,#0
   bNE LOOP                  //   continue loop

   DONE:
   cmp lb,#0              //if(!lb)
   strEQ rb,[curr,gNxt]   //   curr->next=rb;
   cmp rb,#0              //if(!rb)
   strEQ lb,[curr,gNxt]   //   curr->next=lb;
   //epilogue:
   ldmIA SP!,{r4-r6}      //restore caller's variables
bx LR
