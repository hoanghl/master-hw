# $Id: Makefile 12 2008-09-15 20:05:15Z aakurone $ 


CC=gcc -O2 

HEADERS=global.h
OBJECTS= main.o inout.o physical.o neighbourlist.o solve.o forces.o
TARGET= mdmorse


mdmorse:		$(OBJECTS) 
			$(CC) -pg -o $(TARGET) $(OBJECTS) -lm


main.o:			main.c $(HEADERS)
			$(CC) -pg -c main.c

inout.o:		inout.c $(HEADERS)
			$(CC) -pg -c inout.c

physical.o:		physical.c $(HEADERS)
			$(CC) -pg -c physical.c

neighbourlist.o:	neighbourlist.c $(HEADERS)
			$(CC) -pg -c neighbourlist.c

solve.o:		solve.c $(HEADERS)
			$(CC) -pg -c solve.c

forces.o:		forces.c $(HEADERS)
			$(CC) -pg -c forces.c

clean:
			rm -f *.o $(TARGET)

