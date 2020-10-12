#define a0 r0                //holds ptr a0
#define a1 r1                //holds ptr a1
#define diff r2              //holds the addr diff

/*
r0                          r0        r1
unsigned long addrGap(void *a0, void *a1)
{
   if(a0 > a1)
      return a0 - a1;
   return a1 - a0;
}

RETURNS: the positive difference between ptrs r0 and r1
*/

.text
.arm
.globl addrGap
addrGap:
   subs diff,a0,a1           //designate return val to a0-a1
   subLO diff,a1,a0          //if (a0<a1), set return val to a1-a0 instead
   mov r0,diff               //return diff
bx LR  
