CC=gcc
CFLAGS=-marm -Wall -O0 -g

DEP := tstNode.c addrGap.S cntNodes.S findIntMinNode.S

all: $(DEP)
	$(CC) $(CFLAGS) -o xtst $^
