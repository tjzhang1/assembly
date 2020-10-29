#
# Normally, we want to enforce armv4t so that we'll catch it if we use a
# thumb-2 code when we are supposed to be using older ISA.  However, gcc
# doesn't really support this fully anymore, and it causes tools like valgrind
# to freak out, so it's useful to be able to change this temporarily to a
# more modern ISA for debugging, do a make kill, then debug, then put back
# to original value to make sure no thumb-2 code has crept in
#
ARCH=armv4t
CC = gcc
CCFLAGS = -Wall -mthumb -g
CCTFLAGS = -Wall -O3  -march=$(ARCH) -mthumb -mfloat-abi=softfp
# uncomment if needed debug, but don't leave it this way!
#ARCH=armv7-a    # uncomment to debug, but don't leave!
#CCFLAGS = -Wall -mthumb -g
#CCTFLAGS = -Wall -O3  -march=$(ARCH) -mthumb
# Uncomment to compile x[tst,tim]SortC natively for x86
#CC = /usr/bin/arm-linux-gnueabi-gcc
#ARCH=arch
#CCFLAGS = -Wall -g
#CCTFLAGS = -Wall -O3
#
# Use macros to give paths that might change
#
SRCd = ./
#SRCd = /home/whaley/classes/ARMF20/p6_mrgSort
#
# Override this make variable on commandline to try different C/asm codes
#
mrgF=mrgSort

#
# Macros can also help with commonly referenced filenames
#
OBJ = tstSort.o $(mrgF).o
OBJ0 = tstSortC.o $(mrgF).o

force_build:
xtstSortC : tstSortC.o force_build
	$(CC) $(CCFLAGS) -o $@ tstSortC.o $(mrgF).c 
xtimSortC : timSort.o  force_build
	$(CC) $(CCTFLAGS) -o $@ timSort.o $(mrgF).c

xtstSort : $(OBJ)
	$(CC) $(CCFLAGS) -o $@ $(OBJ)
xtimSort : timSort.o mrgSortT.o
	$(CC) $(CCTFLAGS) -o $@ timSort.o mrgSortT.o
#
# This rules remove object and executables when they are no longer needed
#
clean :
	- rm -f $(OBJ) $(OBJ0) timSort.o mrgSortT.o
kill :
	- rm -f $(OBJ) $(OBJ0) timSort.o mrgSortT.o \
                xtstSort xtimSort xtstSortC xtimSortC

tst : xtstSortC xtimSortC
	./xtstSortC 1000000 
	./xtimSortC 1000000
#
# This is an example of a specific rule for building tester.o given tester.c
#
tstSort.o : $(SRCd)/tstSort.c
	$(CC) $(CCFLAGS) -c $(SRCd)/tstSort.c
timSort.o : $(SRCd)/timSort.c
	$(CC) $(CCTFLAGS) -o $@ -c $(SRCd)/timSort.c
#
# Here's an example of renaming the file so we can compile with different
# defines;  $@ means name it whatever is to the left of the : in the rule
#
tstSortC.o : $(SRCd)/tstSort.c
	$(CC) $(CCFLAGS) -DWORD_ONLY -o $@ -c $(SRCd)/tstSort.c

mrgSort.o : mrgSort.S
	$(CC) $(CCFLAGS) -march=$(ARCH) -c mrgSort.S
mrgSortT.o : mrgSort.S
	$(CC) $(CCTFLAGS) -o $@ -march=$(ARCH) -c mrgSort.S
wrdSortT.o : wrdSort.c
	$(CC) $(CCTFLAGS) -march=$(ARCH) -o $@ -c wrdSort.c
#	$(CC) $(CCTFLAGS) -march=armv4t -mfloat-abi=softfp -o $@ -c wrdSort.c
#
# This is a generic rule showing how to build any .o file from a .c file.
# It is only used if no specific rule can be found
#
.c.o:
	$(CC) $(CCFLAGS) -c $<
