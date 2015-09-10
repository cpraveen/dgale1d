# Set the C compiler
CC            = gcc -Wall -O3

HDRS	        = dg.h dg1d.h

LD            = $(CC)

MAKEFILE      = makefile

OBJS          = main.o init.o shape.o gauss.o flux.o atgauss.o \
                update.o div.o result.o project.o time.o \
					 meshvel.o movegrid.o predictor.o

PROGRAM       = dg

all:            $(PROGRAM)

.c.o: $(HDRS)
	$(CC) -c $*.c

$(PROGRAM):     $(OBJS)
				        $(LD) -o $(PROGRAM) $(OBJS) -lm

clean:;         rm -f $(OBJS) core $(PROGRAM)
