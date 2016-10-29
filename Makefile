HEADERS = MSOFT.h MSOFT.h

default: MSOFT.x

MSOFT.o: MSOFT.c $(HEADERS)
	gcc -c MSOFT.c -o MSOFT.o

MSOFT.x: MSOFT.o
	gcc MSOFT.o -o MSOFT.x -lm

clean:
	-rm -f MSOFT.o
	-rm -f MSOFT.x
