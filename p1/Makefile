#
# It's good practice to use make macros for compiler and flags, so when you
# need to change, you change at one place, not throughout makefile
#
CC = gcc
CCFLAGS = -Wall -march=armv7-a -marm -g

#
# Use macros to give paths that might change
#
SRCd = ./
#SRCd = /home/whaley/classes/ARMF20/p1_istr

xitest0 : itester0.o istr.o
	$(CC) $(CCFLAGS) -o xitest0 itester0.o istr.o

xitest : itester.o istr.o
	$(CC) $(CCFLAGS) -o xitest itester.o istr.o


#
# Macros can also help with commonly referenced filenames
#
OBJ = tester.o istr.o

xtest : $(OBJ)
	$(CC) $(CCFLAGS) -o xtest $(OBJ)

#
# This rules remove object and executables when they are no longer needed
#
clean :
	- rm -f $(OBJ) itester.o itester0.o
kill :
	- rm -f $(OBJ) xtest xitest itester.o itester0.o xitest0
#
# This is an example of a specific rule for building tester.o given tester.c
#
tester.o : $(SRCd)/tester.c
	$(CC) $(CCFLAGS) -c $(SRCd)/tester.c
#
# Here's an example of renaming the file so we can compile with different
# defines;  $@ means name it whatever is to the left of the : in the rule
#
itester0.o : $(SRCd)/tester.c
	$(CC) $(CCFLAGS) -DIGNORE_IV=1 -DIGNORE_STR=1 -o $@ -c $(SRCd)/tester.c
itester.o : $(SRCd)/tester.c
	$(CC) $(CCFLAGS) -DIGNORE_STR=1 -o $@ -c $(SRCd)/tester.c

istr.o : istr.S
	$(CC) $(CCFLAGS) -c istr.S
#
# This is a generic rule showing how to build any .o file from a .c file.
# It is only used if no specific rule can be found
#
.c.o:
	$(CC) $(CCFLAGS) -c $<
