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
.thumb
.syntax unified
.local DONE,LOOP
.global mrgChains

mrgChains:
   #define lb r0          //holds the top of left chain
   #define rb r1          //holds the top of right chain
   #define gNxt r2        //holds offset in bytes to ptr to next node
   #define gInt r3        //holds offset in bytes to value of node
   #define vl r4          //holds val of lb (vl=valLeft)
   #define vr r5          //holds val of rb (vr=valRight)
   #define mask r6        //used as a status reg for pseudo-predication
   #define head r7        //holds the head ptr of sorted list
   //prologue:
   push {r4-r7}           //save caller's variables
   movs lb,r0             //save r0 to make room for head
   
   //body:
   ldr vl,[lb,gInt]       //load int value from lb
   ldr vr,[rb,gInt]       //load int value from rb
   
   //initialize head ptr
   /*
   corresponding C code:
      if(valLeft < valRight)
         head = lb;
      else
         head = rb;
   */
   subs mask,vl,vr        //store vl-vr to calculate if (vl-vr)<0
   asrs mask,mask,#31     //mask = (vl-vr < 0) ? 0xFFFFFFFF:0x0
   movs head,lb           //copy lb into head so we don't lose lb's value
   ands head,mask         //head = (vl-vr < 0) ? lb:0x0
   mvns mask,mask         //mask = (vl-vr < 0) ? 0x0:0xFFFFFFFF //(else mask)
   ands mask,rb           //mask = (vl-vr < 0) ? 0x0:rb
   orrs head,mask         //head = (vl-vr < 0) ? lb:rb
   
   mov r12,head           //frees up r7
   #undef head
   #define head r12       //new storage for head ptr is r12
   #define lbnext r7      //r7 now used for loading next value of lb
   //determine which chain (lb or rb) to increment
   /*
   corresponding C code:
      if(valLeft < valRight) {
         lb = lb->next;
         rb = rb;  //must do this explicitly to maintain value
      } else {
         lb = lb;
         rb = rb->next; 
      }
   */
   subs mask,vl,vr        //store vl-vr to calculate if (vl-vr)<0
   asrs mask,mask,#31     //mask = (vl-vr < 0) ? 0xFFFFFFFF:0x0
   ldr lbnext,[lb,gNxt]   //prepare lb->next
   ands lbnext,mask       //ldnext = (vl-vr < 0) ? lb->next:0x0
   ands rb,mask           //rb = (vl-vr < 0) ? rb:0x0

   mvns mask,mask         //mask = (vl-vr < 0) ? 0x0:0xFFFFFFFF //(else mask)
   ands lb,mask           //lb = (vl-vr < 0) ? 0x0:lb
   //make space for rbnext
   #undef vl              //free up r4
   #define rbnext r4      //r4 now used for loading next value of rb
   ldr rbnext,[rb,gNxt]   //prepare rb->next
   ands rbnext,mask       //rbnext = (vl-vr < 0) ? 0x0:rb->next
   orrs lb,lbnext         //lb = (vl-vr < 0) ? lb->next:lb
   orrs rb,rbnext         //rb = (vl-vr < 0) ? rb:rb->next
   /*
   cmp vl,vr              
      ittee LE            //if(valLeft <= valRight)
      movLE head,lb       //   head = lb; (set head of sorted list)
      ldrLE lb,[lb,gNxt]  //   lb=lb->next;
      movGT head,rb       //else head = rb; (set head of sorted list)
      ldrGT rb,[rb,gNxt]  //   rb=rb->next;

   #define curr r7        //holds curr node of sorted list
   movs curr,head          //keep track of current pos w/o losing head ptr

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
   //cleanup; at this pt exactly one chain has nodes, append to end of curr
   DONE:
   cmp lb,#0              
      it EQ               //if(!lb)
      strEQ rb,[curr,gNxt]//   curr->next=rb;
   cmp rb,#0              
      it EQ               //if(!rb)
      strEQ lb,[curr,gNxt]//   curr->next=lb;
   */
   //epilogue:
   pop {r4-r7}            //restore caller's variables
bx LR