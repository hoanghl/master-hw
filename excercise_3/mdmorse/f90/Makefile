# $Id: Makefile 28 2009-07-23 12:30:09Z aakurone $


F90= gfortran -O2

MODULES= modules.o
OBJECTS= main.o inout.o physical.o neighbourlist.o solve.o forces.o
TARGET= mdmorse


mdmorse:		$(OBJECTS) $(MODULES)
			$(F90)  -o $(TARGET) $(OBJECTS) $(MODULES)


main.o:			main.f90 $(MODULES)
			$(F90)  -c main.f90

inout.o:		inout.f90 $(MODULES)
			$(F90)  -c inout.f90

physical.o:		physical.f90 $(MODULES)
			$(F90)  -c physical.f90

neighbourlist.o:	neighbourlist.f90 $(MODULES)
			$(F90)  -c neighbourlist.f90

solve.o:		solve.f90 $(MODULES)
			$(F90)  -c solve.f90

forces.o:		forces.f90 $(MODULES)
			$(F90)  -c forces.f90




modules.o:	modules.f90 
		$(F90)  -c modules.f90


clean:
		rm -f *.o gmon.out *.mod $(TARGET)
