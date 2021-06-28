CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio_1.1.exe : esercizio_1.1.o random.o
	$(CC) random.o esercizio_1.1.o -o esercizio_1.1.exe
esercizio_1.1.o : esercizio_1.1.cpp funzioni.h
	$(CC) -c esercizio_1.1.cpp -o esercizio_1.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	
clean:
	rm *.o esercizio_1.1.exe seed.out
	rm *.dat