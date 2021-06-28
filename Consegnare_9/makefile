CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: esercizio_9.01.exe esercizio_quad_9.01.exe

esercizio_quad_9.01.exe : esercizio_quad_9.01.o random.o
	$(CC) random.o esercizio_quad_9.01.o -o esercizio_quad_9.01.exe
esercizio_9.01.exe : esercizio_9.01.o random.o
	$(CC) random.o esercizio_9.01.o -o esercizio_9.01.exe
esercizio_9.01.o : esercizio_9.01.cpp random.h funzioni.h
	$(CC) -c esercizio_9.01.cpp -o esercizio_9.01.o $(CFLAGS)
esercizio_quad_9.01.o : esercizio_quad_9.01.cpp random.h funzioni.h
	$(CC) -c esercizio_quad_9.01.cpp -o esercizio_quad_9.01.o $(CFLAGS)	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o esercizio_9.01.exe 
	rm *.dat

