CLASSdir = /home/whaley/classes/ARMF20/p7_fparm
#CLASSdir = ./

#CC = gcc
#CC = /usr/bin/arm-linux-gnueabihf-gcc
CCFLAGS = -I$(CLASSdir) -O3 -marm -mfloat-abi=hard -mfpu=vfpv3 -I$(CLASSdir) -g
ASFLAGS = -marm -mfpu=vfpv3 -mfloat-abi=hard -g

objs = dgemmtst.o dgemm_asg.o \
       dgemm_a1b1.o dgemm_a1b0.o dgemm_a1bX.o \
       dgemm_aXb1.o dgemm_aXb0.o dgemm_aXbX.o 

xddottst : ddot.o ddottest.o 
	$(CC) $(CCFLAGS) -o $@ ddottest.o ddot.o 
ddot.o : ddot.S
	$(CC) $(ASFLAGS) -c ddot.S
xsdottst : sdottest.o sdot.o
	$(CC) $(CCFLAGS) -o $@ sdottest.o sdot.o 
#*********************************************
xsdottstC : sdottest.o sdotC.o
	$(CC) $(CCFLAGS) -o $@ sdottest.o sdot.o 
sdotC.o : sdotC.c
	$(CC) $(ASFLAGS) -c sdotC.c
#*********************************************
xsdottstS: sdottest.o sdot_skew.o
	$(CC) $(CCFLAGS) -o $@ sdottest.o sdot_skew.o 
sdot.o : sdot.S
	$(CC) $(ASFLAGS) -c sdot.S
sdot_skew.o : sdot_skew.S
	$(CC) $(ASFLAGS) -c sdot_skew.S

sdottest.o :  $(CLASSdir)/dottest.c
	$(CC) $(CCFLAGS) -DSREAL=1 -c -o $@ $(CLASSdir)/dottest.c
ddottest.o :  $(CLASSdir)/dottest.c
	$(CC) $(CCFLAGS) -DDREAL=1 -c -o $@ $(CLASSdir)/dottest.c
dgemmtst.o :  $(CLASSdir)/gemmtst.c
	$(CC) $(CCFLAGS) -DDREAL=1 -c -o $@ $(CLASSdir)/gemmtst.c

xdgemmtst : $(objs)
	$(CC) $(CCFLAGS) -o $@ $(objs)

dgemm_asg.o : dgemm_asg.S
	$(CC) $(ASFLAGS) -c dgemm_asg.S

dgemm_a1b1.o : $(CLASSdir)/dgemm.c
	$(CC) $(CCFLAGS) -DALPHA1 -DBETA1 -c -o $@ $(CLASSdir)/dgemm.c
dgemm_a1b0.o : $(CLASSdir)/dgemm.c
	$(CC) $(CCFLAGS) -DALPHA1 -DBETA0 -c -o $@ $(CLASSdir)/dgemm.c
dgemm_a1bX.o : $(CLASSdir)/dgemm.c
	$(CC) $(CCFLAGS) -DALPHA1 -DBETAX -c -o $@ $(CLASSdir)/dgemm.c

dgemm_aXb1.o : $(CLASSdir)/dgemm.c
	$(CC) $(CCFLAGS) -DALPHAX -DBETA1 -c -o $@ $(CLASSdir)/dgemm.c
dgemm_aXb0.o : $(CLASSdir)/dgemm.c
	$(CC) $(CCFLAGS) -DALPHAX -DBETA0 -c -o $@ $(CLASSdir)/dgemm.c
dgemm_aXbX.o : $(CLASSdir)/dgemm.c
	$(CC) $(CCFLAGS) -DALPHAX -DBETAX -c -o $@ $(CLASSdir)/dgemm.c
clean :
	rm -f *.o
kill :
	$(MAKE) clean
	rm -f x?dottst xdgemmtst

submit:
	make kill
	$(profd)/bin/p_copy $(PRG)
verify:
	$(profd)/bin/verify $(PRG)
