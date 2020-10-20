#
# It's good practice to use make macros for compiler and flags, so when you
# need to change, you change at one place, not throughout makefile
#
CC = gcc    # for use on a15
#CC = /usr/bin/arm-linux-gnueabihf-gcc
#CC = /usr/bin/arm-linux-gnueabi-gcc  # for use x86 cross-comp ARM
CCF4 = -march=armv4t -mthumb -mfloat-abi=softfp
CCF7 = -mthumb -march=armv7-a
CCFA = -march=armv7-a -marm
CCTFLAGS = -Wall -O3
CCFLAGS = -Wall -g

#
# Use macros to give paths that might change
#
SRCd = ./
#SRCd = /home/whaley/classes/ARMF20/p5_predMerge

CLNlst = tstSort[0,4,A,7].o mrgSort[4,7,A].o mrgChains[4,7,A,F].o \
         timSort[4,7,A].o mrgSort[4,7,A]T.o mrgChains[4,7,A,F]T.o
KILlst = xtest[4,7,A,F] xtime[4,7,A,F] $(CLNlist)
#
# Macros can also help with commonly referenced filenames
#
xtest4 : tstSort4.o mrgSort4.o mrgChains4.o
	$(CC) $(CCFLAGS) $(CCF4) -o $@ tstSort4.o mrgSort4.o mrgChains4.o
xtime4 : timSort4.o mrgSort4T.o mrgChains4T.o
	$(CC) $(CCTFLAGS) $(CCF4) -o $@ timSort4.o mrgSort4T.o mrgChains4T.o
xtest7 : tstSort7.o mrgSort7.o mrgChains7.o
	$(CC) $(CCFLAGS) $(CCF7) -o $@ tstSort7.o mrgSort7.o mrgChains7.o
xtime7 : timSort7.o mrgSort7T.o mrgChains7T.o
	$(CC) $(CCTFLAGS) $(CCF7) -o $@ timSort7.o mrgSort7T.o mrgChains7T.o
xtestA : tstSortA.o mrgSortA.o mrgChainsA.o
	$(CC) $(CCFLAGS) $(CCFA) -o $@ tstSortA.o mrgSortA.o mrgChainsA.o
xtimeA : timSortA.o mrgSortAT.o mrgChainsAT.o
	$(CC) $(CCTFLAGS) $(CCFA) -o $@ timSortA.o mrgSortAT.o mrgChainsAT.o
xtestF : tstSortA.o mrgSortA.o mrgChainsF.o
	$(CC) $(CCFLAGS) $(CCFA) -o $@ tstSortA.o mrgSortA.o mrgChainsF.o
xtimeF : timSortA.o mrgSortAT.o mrgChainsFT.o
	$(CC) $(CCTFLAGS) $(CCFA) -o $@ timSortA.o mrgSortAT.o mrgChainsFT.o

mrgSort4.o : $(SRCd)/mrgSort.c
	$(CC) $(CCFLAGS) $(CCF4) -c -o $@ $(SRCd)/mrgSort.c
mrgSort4T.o : $(SRCd)/mrgSort.c
	$(CC) $(CCTFLAGS) $(CCF4) -c -o $@ $(SRCd)/mrgSort.c
timSort4.o : $(SRCd)/timSort.c
	$(CC) $(CCTFLAGS) $(CCF4) -c -o $@ $(SRCd)/timSort.c
tstSort4.o : $(SRCd)/tstSort.c
	$(CC) $(CCFLAGS) $(CCF4) -c -o $@ $(SRCd)/tstSort.c
mrgChains4.o : mrgChains4.S
	$(CC) $(CCFLAGS) $(CCF4) -c mrgChains4.S
mrgChains4T.o : mrgChains4.S
	$(CC) $(CCTFLAGS) $(CCF4) -o $@ -c mrgChains4.S

mrgSort7.o : $(SRCd)/mrgSort.c
	$(CC) $(CCFLAGS) $(CCF7) -c -o $@ $(SRCd)/mrgSort.c
mrgSort7T.o : $(SRCd)/mrgSort.c
	$(CC) $(CCTFLAGS) $(CCF7) -c -o $@ $(SRCd)/mrgSort.c
timSort7.o : $(SRCd)/timSort.c
	$(CC) $(CCTFLAGS) $(CCF7) -c -o $@ $(SRCd)/timSort.c
tstSort7.o : $(SRCd)/tstSort.c
	$(CC) $(CCFLAGS) $(CCF7) -c -o $@ $(SRCd)/tstSort.c
mrgChains7.o : mrgChains7.S
	$(CC) $(CCFLAGS) $(CCF7) -c mrgChains7.S
mrgChains7T.o : mrgChains7.S
	$(CC) $(CCTFLAGS) $(CCF7) -o $@ -c mrgChains7.S

mrgSortA.o : $(SRCd)/mrgSort.c
	$(CC) $(CCFLAGS) $(CCFA) -c -o $@ $(SRCd)/mrgSort.c
mrgSortAT.o : $(SRCd)/mrgSort.c
	$(CC) $(CCTFLAGS) $(CCFA) -c -o $@ $(SRCd)/mrgSort.c
timSortA.o : $(SRCd)/timSort.c
	$(CC) $(CCTFLAGS) $(CCFA) -c -o $@ $(SRCd)/timSort.c
tstSortA.o : $(SRCd)/tstSort.c
	$(CC) $(CCFLAGS) $(CCFA) -c -o $@ $(SRCd)/tstSort.c
mrgChainsA.o : mrgChainsA.S
	$(CC) $(CCFLAGS) $(CCFA) -c mrgChainsA.S
mrgChainsAT.o : mrgChainsA.S
	$(CC) $(CCTFLAGS) $(CCFA) -o $@ -c mrgChainsA.S
mrgChainsF.o : mrgChainsF.S
	$(CC) $(CCFLAGS) $(CCFA) -c mrgChainsF.S
mrgChainsFT.o : mrgChainsF.S
	$(CC) $(CCTFLAGS) $(CCFA) -o $@ -c mrgChainsF.S

xtime : timSort.o mrgSortT.o
	$(CC) $(CCTFLAGS) -o $@ timSort.o mrgSortT.o
#
# This rules remove object and executables when they are no longer needed
#
clean:
	rm -f $(CLNlst)
kill :
	rm -f $(CLNlst) $(KILlst)
#
# Here's an example of renaming the file so we can compile with different
# defines;  $@ means name it whatever is to the left of the : in the rule
#
tstSort0.o : $(SRCd)/tstSort.c
	gcc $(CCFLAGS) -DWORD_ONLY -o $@ -c $(SRCd)/tstSort.c

# This is a generic rule showing how to build any .o file from a .c file.
# It is only used if no specific rule can be found
#
#.c.o:
#	$(CC) $(CCFLAGS) -c $<
