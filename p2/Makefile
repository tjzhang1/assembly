#
# It's good practice to use make macros for compiler and flags, so when you
# need to change, you change at one place, not throughout makefile
#
CC = gcc
CCFLAGS = -Wall -marm -g

#
# Use macros to give paths that might change
#
SRCd = ./
asm=lvecsumL
#SRCd = /home/whaley/classes/ARMF20/p2_isum
#
# Macros can also help with commonly referenced filenames
#
OBJ = itester.o $(asm).o


force_build :
xitest : force_build $(OBJ)
	$(CC) $(CCFLAGS) -o xitest $(OBJ)


#
# This rules remove object and executables when they are no longer needed
#
clean :
	- rm -f $(OBJ)
kill :
	- rm -f $(OBJ) xitest
#
# Here's an example of renaming the file so we can compile with different
# defines;  $@ means name it whatever is to the left of the : in the rule
#
itester.o : $(SRCd)/lsum_test.c
	$(CC) $(CCFLAGS) -o $@ -c $(SRCd)/lsum_test.c

$(asm).o : $(asm).S
	$(CC) $(CCFLAGS) -c $(asm).S
#
# This is a generic rule showing how to build any .o file from a .c file.
# It is only used if no specific rule can be found
#
.c.o:
	$(CC) $(CCFLAGS) -c $<
