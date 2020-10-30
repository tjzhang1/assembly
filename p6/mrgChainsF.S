/*
r0
void *mrgChains
(
   void *lb,              //r0
   void *rb,              //r1
   uint gNxt,             //r2
   uint gInt,             //r3
   void *endL,            //SP[0]
   void *endR,            //SP[4]
   void *ML               //SP[8]
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
.syntax unified
.local DONE1,DONE2,LOOP,ELSE1,ELSE2,QRET1,QRET2
.global mrgChains

#define FSZ 4*6        //6 regs stored
#define arg4 FSZ       //offset of first arg on SP
#define arg5 FSZ+4     //offset of 2nd arg on SP
#define arg6 FSZ+8     //offset of 3rd arg on SP

mrgChains:
   #define lb r0          //holds the top of left chain
   #define rb r1          //holds the top of right chain
   #define gNxt r2        //holds offset in bytes to ptr to next node
   #define gInt r3        //holds offset in bytes to value of node
   #define vl r4          //holds first val of lb = min val of lb
   #define maxR r5        //holds val of endR = max val of rb
   #define endL r6        //holds last node of lb
   #define endR r7        //holds last node of rb
   push {r4-r9}           //save caller's variables
   ldr endL,[SP,arg4]     //load endL argument from stack
   ldr endR,[SP,arg5]     //load endR argument from stack
   /*
      if(maxL <= rbVal)
      {
         *ML = endR;
         endL->next=rb;
         return lb;
      }
      else if (maxR <= lbVal)
      {
         *ML = endL;
         endR->next=lb;
         return rb;
      }
   */
   ldr vl,[lb,gInt]       //load int value from lb
   ldr maxR,[endR,gInt]   //load int value from endR
   cmp maxR,vl            //if(max(rb) <= min(lb))
   bLE QRET1              //   concatenate lb to rb and return
   mov r8,vl              //save vl by moving to higher reg
   #undef vl              //reassign vl for vr
   #undef maxR            //free r5 for maxL
   #define vr r4          //holds first val of rb = min val of rb
   #define maxL r5        //holds val of endL = max val of lb
   ldr vr,[rb,gInt]       //load int value from rb
   ldr maxL,[endL,gInt]   //load int value from endL
   cmp maxL,vr            //if(max(lb) <= min(rb))
   bLE QRET2              //   concatenate rb to lb and return
   mov r5,vl              //restore vl from high reg
   mov r8,endL            //store endL to high reg
   mov r9,endR            //store endR to high reg
   mov r7,lb              //make room for new list head
   #undef maxL            //free r5 for vl
   #undef endL            //free r6 for curr
   #undef endR            //free r7 for lb
   #undef lb              //reassign lb
   #define head r0        //holds the base ptr of sorted list
   #define vl r5          //new spot for vl
   #define curr r6        //holds curr node of sorted list
   #define lb r7          //new spot for lb
   #define endL r8        //new spot for endL
   #define endR r9        //new spot for endR
   //new reg configuration
   //r0: head
   //r1: rb
   //r2: gNxt
   //r3: gInt
   //r4: vr
   //r5: vl
   //r6: curr
   //r7: lb
   //r8: endL
   //r9: endR
   /*
      Initialize base pointer head
   */
   cmp vl,vr              //if(valLeft <= valRight)
      bGT ELSE1
      movs head,lb        //   head = lb; (set head of sorted list)
      ldr lb,[lb,gNxt]    //   lb=lb->next; (increment lb)
      movs curr,head      //   keep track of current pos w/o losing head ptr
      cmp lb,#0           //   if(!lb)
      bEQ DONE1           //      skip loop
      ldr vl,[lb,gInt]    //   else update int value from lb
      b LOOP              //   begin loop
   ELSE1:                 //else if(valLeft > valRight)
      movs head,rb        //   head = rb; (set head of sorted list)
      ldr rb,[rb,gNxt]    //   rb=rb->next; (increment rb)
      movs curr,head      //   keep track of current pos w/o losing head ptr
      cmp rb,#0           //   if(!rb)
      bEQ DONE2           //      skip loop
      ldr vr,[rb,gInt]    //   else update int value from rb, begin loop

   LOOP:
      cmp vl,vr              //if(valLeft <= valRight)
         bGT ELSE2
         str lb,[curr,gNxt]  //   set nxt node in sorted list = lb
         movs curr,lb        //   curr = lb; (curr=curr->next)
         ldr lb,[lb,gNxt]    //   lb=lb->next;
         cmp lb,#0           //   if(!lb)
         bEQ DONE1           //      break;
         ldr vl,[lb,gInt]    //   else update int value from lb
         b LOOP              //      continue;
      ELSE2:                 //else if(valLeft > valRight)
         str rb,[curr,gNxt]  //   set nxt node in sorted list = rb
         movs curr,rb        //   curr = rb; (curr=curr->next)
         ldr rb,[rb,gNxt]    //   rb=rb->next;
         cmp rb,#0           //   if(!rb)
         bEQ DONE2           //      break;
         ldr vr,[rb,gInt]    //   else update int value from rb
   b LOOP                    //      continue;


   //epilogue:
   pop {r4-r9}      //restore caller's variables
bx LR


QRET1:  //occurs if(max(rb) <= min(lb))
   #define lb r0          //holds the top of left chain
   #define rb r1          //holds the top of right chain
   #define gNxt r2        //holds offset in bytes to ptr to next node
   #define ML r4          //updated to point to last node
   #define endL r6        //holds last node of lb
   #define endR r7        //holds last node of rb
   ldr ML,[SP,arg6]       //load ML from stack
   str endL,[ML]          //*ML = endL
   str lb,[endR,gNxt]     //append lb to end of rb
   movs r0,rb             //return rb
   pop {r4-r9}      //restore caller's variables
bx LR

QRET2:  //occurs if(max(lb) <= min(rb))
   #define lb r0          //holds the top of left chain
   #define rb r1          //holds the top of right chain
   #define gNxt r2        //holds offset in bytes to ptr to next node
   #define ML r4          //updated to point to last node
   #define endL r6        //holds last node of lb
   #define endR r7        //holds last node of rb
   ldr ML,[SP,arg6]       //load ML from stack
   str endR,[ML]          //*ML = endL
   str rb,[endL,gNxt]     //append rb to end of lb
   pop {r4-r9}      //restore caller's variables
bx LR                     //return lb

DONE1:  //occurs if(!lb)
   #define rb r1          //holds the top of right chain
   #define gNxt r2        //holds offset in bytes to ptr to next node
   #define ML r4          //updated to point to last node
   #define curr r6        //holds curr node of sorted list
   #define endR r9        //holds last node of rb
   ldr ML,[SP,arg6]       //load ML from stack
   mov r5,endR            //mov endR from higher register
   str r5,[ML]            //*ML = endR
   str rb,[curr,gNxt]     //append rb to current node
bx LR

DONE2:  //occurs if(!rb)
   #define gNxt r2        //holds offset in bytes to ptr to next node
   #define ML r4          //updated to point to last node
   #define curr r6        //holds curr node of sorted list
   #define lb r7          //holds the top of left chain
   #define endL r8        //holds last node of lb
   ldr ML,[SP,arg6]       //load ML from stack
   mov r5,endL            //mov endL from higher register
   str r5,[ML]            //*ML = endL
   str lb,[curr,gNxt]     //append lb to current node
bx LR
