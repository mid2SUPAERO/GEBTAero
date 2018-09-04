include gebtAeroMake # Change this for specific program

OBJS     = $(OBJS90)

COMPILE  = gfortran -I/usr/include -I/usr/include/mumps_seq -static -static-libgfortran -static-libgcc 
LINK     = gfortran -static-libgfortran
cFLAGS   = -O3 -fbounds-check -c 
lFLAGS   = -O3 -fbounds-check
LIBS     = -llapack -larpack -lmpiseq_seq -ldmumps_seq
SRCDIR   =../src/

$(PGM):  $(OBJS)   
	$(LINK) $(lFLAGS) -o $(PGM)   $(OBJS) $(LIBS)
	rm -rf *.o *.c *.mod

$(OBJS90): %.o: $(SRCDIR)%.f90
	$(COMPILE) $(cFLAGS) $<

$(OBJSF): %.o: $(SRCDIR)%.f 
	$(COMPILE) $(cFLAGS) $<

.PHONY: clean,unical
clean :
	rm -rf *.o *.c *.mod gebtaero unical
	
unical :
	gcc -o unical.o -c ../src/unical/unical1_b.c 
	gcc -o unical unical.o
	rm unical.o
