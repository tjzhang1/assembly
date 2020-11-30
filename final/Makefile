CC = gcc
#uncomment to enable ARM assembly cross-compilation on x86 (if installed)
#CC = /usr/bin/arm-linux-gnueabi-gcc
#
ARCH=armv7-a
CCFLAGS = -Wall -march=$(ARCH) -mthumb -g
CCTFLAGS = -Wall -O3  -march=$(ARCH) -mthumb
# Uncomment to compile x[tst,tim]SortC natively for x86
#ARCH=arch
#CCFLAGS = -Wall -g
#CCTFLAGS = -Wall -O3
#
# Use macros to give paths that might change
#
SRCd = ./
#SRCd = /home/whaley/classes/ARMF20/pF_mrgSortUL
#
# Override this make variable on commandline to try different C/asm codes
#
mrgF=mrgSort

#
# Macros can also help with commonly referenced filenames
#

VG = valgrind --leak-check=full --malloc-fill=0x88 --track-origins=yes
#uncomment if on x86 wt emulation
#VG = qemu-arm -L /usr/arm-linux-gnueabihf
# sizes to run test for
SZS = 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 \
      26 27 777 5001

TSZS = 1000000 1000000
testC: xtstSortC 
	$(VG) ./xtstSortC $(SZS)
testC0: xtstSortC0
	$(VG) ./xtstSortC0 $(SZS)
timeC: xtimSortC 
	./xtimSortC $(TSZS)
timeC0: xtimSortC0 
	./xtimSortC0 $(TSZS)
test : xtstSort
	$(VG) ./xtstSort $(SZS)
test0: xtstSort0
	$(VG) ./xtstSort0 $(SZS)
time: xtimSort 
	./xtimSort $(TSZS)
time0: xtimSort0 
	./xtimSort0 $(TSZS)

force_build:
xtstSortC : force_build
	$(CC) $(CCFLAGS) -DWORD_ONLY -o $@ tstSort.c $(mrgF).c
xtimSortC : force_build
	$(CC) $(CCTFLAGS) -DWORD_ONLY -o $@ timSort.c $(mrgF).c
xtstSortC0 : force_build
	$(CC) $(CCFLAGS) -DWORD_ONLY -DSMALL_WORD -o $@ tstSort.c $(mrgF).c
xtimSortC0 : force_build
	$(CC) $(CCTFLAGS) -DWORD_ONLY -DSMALL_WORD -o $@ tstSort.c $(mrgF).c

xtstSort : force_build
	$(CC) $(CCFLAGS) -o $@ tstSort.c $(mrgF).S
xtimSort :  force_build
	$(CC) $(CCTFLAGS) -o $@ timSort.c $(mrgF).S
xtstSort0 : force_build
	$(CC) $(CCFLAGS) -o $@ tstSort.c $(mrgF).S
xtimSort0 :  force_build
	$(CC) $(CCTFLAGS) -o $@ timSort.c $(mrgF).S
#
# This rules remove object and executables when they are no longer needed
#
clean :
	- rm -f $(OBJ) $(OBJ0) timSort.o mrgSortT.o
kill :
	- rm -f $(OBJ) $(OBJ0) timSort.o mrgSortT.o \
                xtstSort xtimSort xtstSortC xtimSortC
#
# This is a generic rule showing how to build any .o file from a .c file.
# It is only used if no specific rule can be found
#
.c.o:
	$(CC) $(CCFLAGS) -c $<
