#
# Override this macro to try a different sort
#
mySort=selSort
#
# It's good practice to use make macros for compiler and flags, so when you
# need to change, you change at one place, not throughout makefile
#
#CC = /usr/bin/arm-linux-gnueabi-gcc
#CCFLAGS = -Wall -mthumb -g
#CCTFLAGS = -Wall -O3  -march=armv4t -mthumb -mfloat-abi=softfp
CC = gcc
CCFLAGS = -Wall -g
CCTFLAGS = -Wall -O3 

#
# Use macros to give paths that might change
#
SRCd = ./
#SRCd = ~rcwhaley/classFiles/ARMF20/p0_mrgSort

xtstSort : force_build
	rm -f mySortG.o tstSort.o
	$(MAKE) mySortG.o tstSort.o mySort=$(mySort)
	$(CC) $(CCFLAGS) -o $@ tstSort.o mySortG.o

xtimSort : force_build
	rm -f mySortT.o timSort.o
	$(MAKE) mySortT.o timSort.o mySort=$(mySort)
	$(CC) $(CCTFLAGS) -o $@ timSort.o mySortT.o

#
# Target that doesn't exist, with no build will force rebuild
#
force_build:
#
# This rules remove object and executables when they are no longer needed
#
clean :
	- rm -f tstSort.o timSort.o mySortT.o mySortG.o
kill :
	make clean
	rm -f xtstSort xtimSort
#
# This is an example of a specific rule for building tester.o given tester.c
#
tstSort.o : $(SRCd)/tstSort.c
	$(CC) $(CCFLAGS) -DmySort=$(mySort) -o $@ -c $(SRCd)/tstSort.c
timSort.o : $(SRCd)/timSort.c
	$(CC) $(CCFLAGS) -DmySort=$(mySort) -o $@ -c $(SRCd)/timSort.c
#
# Here's an example of renaming the file so we can compile with different
# defines;  $@ means name it whatever is to the left of the : in the rule
#
mySortT.o :
	$(CC) $(CCTFLAGS) -DmySort=$(mySort) -o $@ -c $(mySort).c
mySortG.o :
	$(CC) $(CCFLAGS) -DmySort=$(mySort) -o $@ -c $(mySort).c
#
# This is a generic rule showing how to build any .o file from a .c file.
# It is only used if no specific rule can be found
#
.c.o:
	$(CC) $(CCFLAGS) -c $<

slow:
	make xtstSort mySort=mrgSortSlow
	make xtimSort mySort=mrgSortSlow

fast:
	make xtstSort mySort=mrgSortFast
	make xtimSort mySort=mrgSortFast

test:
	valgrind --leak-check=full --malloc-fill=0x88 --track-origins=yes \./xtstSort 0 1 2 3 4 5 6 7 8 9 10 11 777

tim:
	./xtimSort 5000000 5000000 5000000
