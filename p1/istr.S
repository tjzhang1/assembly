#define i1 r0        //operand for sum
#define i2 r1        //operand sum
#define i3 r2        //operand for sum and multiply
#define pi4 r3       //ptr to operand for sum and multiply
#define i4 r12       //holds value from pi4
#define N5 r4        //holds #5 for multiply
#define indx r5      //holds index for for loop
#define ALP r6       //holds alp[n]
#define pALP r7      //holds ptr char *alp
#define FSZ 4*4      //saving caller's r4,r5,r6,r7
/*
 *             r0      r1      r2       r3        SP
 *int istr(int i1, int i2, int i3, int *i4, char *alp)
 *{
 *   char n;
 *   for(n=0; n<26; n++)
 *      alp[n] = (char) (n+97);
 *   alp[n] = 0;
 *   *i4 = 5 * (*i4) + i3;
 *   return i1 + i2;
 *}
 */
.text
.arm
.local LOOP1
.global istr
istr: //prologue
   sub SP, SP, #FSZ      //allocate stack space
   str r4, [SP]          //store caller's r4
   str r5, [SP,#4]       //store caller's r5
   str r6, [SP,#8]       //store caller's r6
   str r7, [SP,#12]      //store caller's r7
   ldr pALP,[SP,#FSZ]    //load char *alp (arg 5) from stack
/*
 *for(n=0; n<26; n++)
 *   alp[n] = (char) (n+97);
 *alp[n] = 0;
 */   
   add pALP,pALP,#26     //pALP += 26
   mov indx, #-26        //n = -26
LOOP1:
   add ALP,indx,#123     //ALP = n+123: n=[-26,0)-->ALP=[97,123)
   strb ALP,[pALP,indx]  //pALP[n] = ALP
   adds indx,indx,#1     //n++
bNE LOOP1                //go to LOOP1
   mov ALP, #0           //ALP = '\0'
   strb ALP,[pALP,indx]  //pALP[26] = ALP

   ldr i4, [pi4]         //ld value pointed to by pi4
   mov N5, #5            //store 5 into reg for multiply
   mla i4,N5,i4,i3       //*i4 = 5*(*i4) + i3
   str i4, [pi4]         //point pi4 to newly calculated value

   add i1, i1, i2        //return i1 + i2

//epilogue - return and restore regs and frame
   ldr r4, [SP]          //restore caller's r4
   ldr r5, [SP,#4]       //restore caller's r5
   ldr r6, [SP,#8]       //restore caller's r6
   ldr r7, [SP,#12]      //restore caller's r7
   add SP, SP, #FSZ      //dealloc stack space
   bx lr
