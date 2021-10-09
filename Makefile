#Makefile
CC	= gcc
CFLAGS	= -Wall -O2 -I/usr/local/include
LIBS = -lm
OBJS = main.o unit.o run.o calc_energy.o gaussian.o setup.o
PROGRAM = nbody.out

all: $(PROGRAM)
clean:; rm $(OBJS) nbody.out 

$(PROGRAM):	$(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(PROGRAM)

main.o:	main.c
	$(CC) -c main.c
unit.o: unit.c
	$(CC) -c unit.c
run.o: run.c
	$(CC) -c run.c
calc_energy.o: calc_energy.c
	$(CC) -c calc_energy.c
gaussian.o: gaussian.c
	$(CC) -c gaussian.c
setup.o: setup.c
	$(CC) -c setup.c


